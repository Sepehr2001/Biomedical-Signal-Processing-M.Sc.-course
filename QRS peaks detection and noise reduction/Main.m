%% ECG Noise Reduction using Synchronous Averaging 
clear; close all; clc;

%% ===== STEP 0: Load ECG Data =====
ecg_raw = load('ecg.dat');
fs = 1000;           % Sampling frequency (Hz)
N = length(ecg_raw);
t = (0:N-1) / fs;    % Time vector (seconds)

%% ===== STEP 1: Visualize Raw Noisy Signal =====
figure('Name','Step 1: Raw Noisy ECG');
plot(t, ecg_raw);           xlabel('Time (s)');             ylabel('Amplitude (mV)');
title('Step 1: Raw Noisy ECG Signal', 'FontSize', 14, 'FontWeight', 'bold');        grid on;

% %% ===== STEP 2: R-Peak Detection (Pan-Tompkins based on built in Matlab functions ) =====
% 
% % Bandpass filter: 5-15 Hz to enhance QRS
% [b_bp, a_bp] = butter(3, [5 15]/(fs/2), 'bandpass');
% ecg_filtered = filtfilt(b_bp, a_bp, ecg_raw);
% 
% % Derivative
% ecg_diff = diff(ecg_filtered);
% ecg_diff = [ecg_diff(1); ecg_diff];
% 
% % Squaring
% ecg_sq = ecg_diff .^ 2;
% 
% % Moving window integration (150 ms window)
% win_size = round(0.15 * fs);
% ecg_mwi = movmean(ecg_sq, win_size);
% 
% % Normalize
% ecg_mwi = ecg_mwi / max(ecg_mwi);
% 
% % Threshold for R-peak detection
% threshold = 0.35;
% min_rr_samples = round(0.4 * fs);   % Minimum RR: 400ms (150 bpm max)  {60 s/150 bpm} = 0.4s
% 
% [~, R_locs] = findpeaks(ecg_mwi, 'MinPeakHeight', threshold, ...
%                                  'MinPeakDistance', min_rr_samples);
% 
% % Remove duplicate R-peaks
% R_locs = unique(R_locs);
% num_beats = length(R_locs);
% clc
% fprintf('Number of detected R-peaks: %d\n', num_beats);
% fprintf('Average Heart Rate: %.1f bpm\n', (num_beats / t(end)) * 60);
% 
% % ---- Plot R-peak detection pipeline ----
% figure('Name','Step 2: R-Peak Detection');
% subplot(4,1,1);
% plot(t, ecg_raw, 'Color',[0.2 0.4 0.8]);
% title('ًRaw ECG + Detected R-peaks','FontSize',11,'FontWeight','bold'); 
% ylabel('Amplitude'); grid on;
% hold on;
% plot(t(R_locs), ecg_raw(R_locs), 'rv', 'MarkerFaceColor','r','MarkerSize',8);
% 
% subplot(4,1,2);
% plot(t, ecg_filtered, 'Color',[0.1 0.6 0.3]);
% title('Bandpass Filtered (5-15 Hz)','FontSize',11,'FontWeight','bold'); ylabel('Amplitude'); grid on;
% 
% subplot(4,1,3);
% plot(t, ecg_sq, 'Color',[0.8 0.4 0.1]);
% title('Squared Derivative (Energy)','FontSize',11,'FontWeight','bold'); ylabel('Energy'); grid on;
% 
% subplot(4,1,4);
% plot(t, ecg_mwi, 'Color',[0.6 0.1 0.6],'LineWidth',1.2); hold on;
% yline(threshold, 'r--', 'LineWidth',1.5, 'Label','Threshold');
% plot(t(R_locs), ecg_mwi(R_locs), 'rv', 'MarkerFaceColor','r','MarkerSize',8);
% title('Moving Window Integration + Detected R-peaks','FontSize',11,'FontWeight','bold');
% ylabel('Normalized'); xlabel('Time (s)'); grid on;
% for k = 1:4; subplot(4,1,k); xlim([t(1) t(end)]); end

%% ===== STEP 2.2: R-Peak Detection (Pan-Tompkins based on Dr. Moradi's slides and Rangayyan Book) =====
clc
out_pan_tomp = pan_tomp(ecg_raw, fs, 1);

% Threshold for R-peak detection
threshold = 0.35;
min_rr_samples = round(0.4 * fs);   % Minimum RR: 400ms (150 bpm max)  {60 s/150 bpm} = 0.4s

[~, R_locs] = findpeaks(out_pan_tomp, 'MinPeakHeight', threshold, ...
                                 'MinPeakDistance', min_rr_samples);

% Remove duplicate R-peaks
R_locs = unique(R_locs);
num_beats = length(R_locs);
clc
fprintf('Number of detected R-peaks: %d\n', num_beats);
fprintf('Average Heart Rate: %.1f bpm\n', (num_beats / t(end)) * 60);

% ---- Plot R-peak detection pipeline ----
figure('Name','Step 2: R-Peak Detection');
subplot(2,1,1);
plot(t, ecg_raw, 'Color',[0.2 0.4 0.8]);
title('ًRaw ECG + Detected R-peaks','FontSize',11,'FontWeight','bold'); 
ylabel('Amplitude'); grid on;
hold on;
plot(t(R_locs), ecg_raw(R_locs), 'rv', 'MarkerFaceColor','r','MarkerSize',8);

subplot(2,1,2);
plot(t, out_pan_tomp, 'Color',[0.6 0.1 0.6],'LineWidth',1.2); hold on;
yline(threshold, 'r--', 'LineWidth',1.5, 'Label','Threshold');
plot(t(R_locs), out_pan_tomp(R_locs), 'rv', 'MarkerFaceColor','r','MarkerSize',8);
title('Moving Window Integration + Detected R-peaks','FontSize',11,'FontWeight','bold');
ylabel('Normalized'); xlabel('Time (s)'); grid on;
for k = 1:2; subplot(2,1,k); xlim([t(1) t(end)]); end



%% ===== STEP 2: Beat Segmentation (R-to-R based) =====
clc

% RR intervals (in samples)
RR_intervals = diff(R_locs);   % length = num_beats - 1. this is the R-R lengths

% Each beat: from current R-peak to next R-peak
% So we get num_beats-1 beats (last beat has no "next R")
num_valid    = length(R_locs) - 1;
R_locs_valid = R_locs(1:num_valid);

fprintf('Valid beats: %d\n', num_valid);
fprintf('RR intervals (samples): min=%d, max=%d, mean=%.1f\n', ...
        min(RR_intervals), max(RR_intervals), mean(RR_intervals));

% Store beats as cell array (different lengths)
beats_cell = cell(1, num_valid);
for i = 1:num_valid
    win           = R_locs(i) : R_locs(i+1) - 1;   % R_k to R_{k+1} (exclusive)
    beats_cell{i} = ecg_raw(win);
end

% ── Plot all beats overlaid 
figure('Name', 'Step 3: Beat Segmentation (R-to-R)');
hold on;
colors = lines(num_valid);
for i = 1:num_valid
    beat   = beats_cell{i};
    plot( beat, 'Color', [colors(i,:) 0.5]);
end

xlabel('samples');
ylabel('Amplitude');
title(sprintf('Step 2: %d Beats Segmented by R-to-R Distance', num_valid));
grid on;

% ── Print beat lengths ────────────────────────────────────────────────────
fprintf('\nBeat lengths (samples):\n');
for i = 1:num_valid
    fprintf('  Beat %2d: %d samples (%.1f ms)\n', i, length(beats_cell{i}), length(beats_cell{i})/fs*1000);
end

%% ===== STEP 3B: Align beats to median-length beat =====
beat_lengths = cellfun(@length, beats_cell);
[~, ref_idx] = min(abs(beat_lengths - median(beat_lengths)));
ref_len      = beat_lengths(ref_idx);
fprintf('Reference beat: #%d | length: %d samples (%.1f ms)\n', ref_idx, ref_len, ref_len/fs*1000);

ref_beat = beats_cell{ref_idx};
ref_amp  = max(abs(ref_beat));

beats_aligned = zeros(ref_len, num_valid);
scale_len     = zeros(1, num_valid);
scale_amp     = zeros(1, num_valid);

for i = 1:num_valid
    original_beat   = beats_cell{i};
    scale_len(i)    = length(original_beat) / ref_len;
    beat_resampled  = interp1(linspace(0,1,length(original_beat)), original_beat, linspace(0,1,ref_len), 'spline');
    curr_amp        = max(abs(beat_resampled));
    scale_amp(i)    = curr_amp / ref_amp;
    beats_aligned(:,i) = beat_resampled / scale_amp(i);
end

fprintf('\nScale factors (length | amplitude):\n');
for i = 1:num_valid
    fprintf('  Beat %2d: len_scale=%.4f | amp_scale=%.4f\n', i, scale_len(i), scale_amp(i));
end

t_ref = (0:ref_len-1) / fs * 1000;

figure('Name', 'Step 3B: Aligned Beats');
hold on;
for i = 1:num_valid
    if i == ref_idx
        continue
    end
    plot(t_ref, beats_aligned(:,i), 'Color',[0.5 0.7 1 0.5], 'LineWidth', 0.6, 'HandleVisibility','off');
end
h_ref  = plot(t_ref, beats_aligned(:,ref_idx), 'r-',  'LineWidth', 2);
h_mean = plot(t_ref, mean(beats_aligned, 2),    'k--', 'LineWidth', 2);
xlabel('Time (ms)'); ylabel('Amplitude');
title(sprintf('Beats Aligned to Reference Beat #%d (%d samples)', ref_idx, ref_len));
legend([h_ref, h_mean], {'Reference Beat', 'Mean of All Beats'}, 'Location', 'northeast');
grid on; box on;

%% ===== STEP 4 & 5: Loop over num_beats = 2 to 11 =====
max_beats  = min(11, num_valid);
beat_range = 2:max_beats;
num_iter   = 1000;

ref_len      = beat_lengths(ref_idx);
T_samp   = 1/fs;
y_bar    = mean(beats_aligned, 2);
sigma2_y = (1/(ref_len * T_samp)) * sum(y_bar .^ 2);

SNR_dB_matrix = zeros(num_iter, length(beat_range));   % 100 x 10
ecg_all       = zeros(N, length(beat_range));           % store last iteration

for iter = 1:num_iter
    for b_idx = 1:length(beat_range)
        num_beats = beat_range(b_idx);

        % ── Template ──────────────────────────────────────────────────────
        chosen_global    = randperm(num_valid, num_beats);
        template_aligned = mean(beats_aligned(:, chosen_global), 2);

        % ── SNR ───────────────────────────────────────────────────────────
        sigma2_eta = 0;
        for k = 1:num_beats
            for n = 1:ref_len
                sigma2_eta = sigma2_eta + (beats_aligned(n, chosen_global(k)) - y_bar(n))^2;
            end
        end
        sigma2_eta = sigma2_eta / (ref_len * T_samp * (num_beats - 1));
        SNR_dB_matrix(iter, b_idx) = 10 * log10(sigma2_y / sigma2_eta);

        % ── Reconstruct (only save last iteration) ────────────────────────
        if iter == num_iter
            beats_recovered = cell(1, num_valid);
            for i = 1:num_valid
                beats_recovered{i} = interp1(linspace(0,1,ref_len), ...
                                              template_aligned * scale_amp(i), ...
                                              linspace(0,1,beat_lengths(i)), 'spline');
            end
            ecg_reconstructed = zeros(N, 1);
            weight_map        = zeros(N, 1);
            for i = 1:num_valid
                win = R_locs(i) : R_locs(i+1)-1;
                if length(win) == length(beats_recovered{i})
                    ecg_reconstructed(win) = ecg_reconstructed(win) + beats_recovered{i}';
                    weight_map(win)        = weight_map(win) + 1;
                end
            end
            ecg_all(:, b_idx) = ecg_reconstructed;
        end
    end
    fprintf('Iteration %3d / %d done\n', iter, num_iter);
end

% ── Statistics ────────────────────────────────────────────────────────────
SNR_mean = mean(SNR_dB_matrix, 1);   % [1 x 10]
SNR_std  = std(SNR_dB_matrix,  0, 1);

fprintf('\n====== SNR Statistics (100 iterations) ======\n');
fprintf('%-12s %-12s %-10s\n', 'num_beats', 'Mean SNR', 'Std SNR');
fprintf('%s\n', repmat('-',1,36));
for b_idx = 1:length(beat_range)
    fprintf('  %-10d %-12.2f %.2f dB\n', beat_range(b_idx), SNR_mean(b_idx), SNR_std(b_idx));
end

% ── SNR Plot with error bars ───────────────────────────────────────────────
figure('Name', 'SNR vs num_beats (100 iterations)');
errorbar(beat_range, SNR_mean, SNR_std, 'b-o', ...
         'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 7, ...
         'CapSize', 8);
xlabel('Number of Beats Used for Averaging', 'FontSize', 12);
ylabel('SNR (dB)',                           'FontSize', 12);
title('Mean SNR ± Std over 100 Iterations',  'FontSize', 13, 'FontWeight', 'bold');
xticks(beat_range); grid on; box on;

%% ===== Plot: noisy vs reconstructed signals ======
plot_beats   = [2, 5, 11];
plot_indices = plot_beats - 1;   % beat_range starts at 2, so index = num_beats - 1

figure('Name', 'Last Iteration: Reconstructed ECG (2, 5, 11 beats)');
for p = 1:3
    subplot(3, 1, p); hold on;
    b_idx = plot_indices(p);
    plot(t, ecg_raw,           'Color',[0.8 0.3 0.3 0.6], 'LineWidth', 0.6);
    plot(t, ecg_all(:, b_idx), 'Color',[0.1 0.6 0.2],     'LineWidth', 1.2);
    title(sprintf('num\\_beats=%d  |  SNR=%.1f dB', ...
                  plot_beats(p), SNR_dB_matrix(end, b_idx)), 'FontSize', 10);
    xlabel('Time (s)'); ylabel('Amplitude');
    xlim([t(1) t(end)]); grid on; box on;
    if p == 1
        legend('Noisy ECG', 'Reconstructed ECG', 'Location', 'northeast', 'FontSize', 8);
    end
end
sgtitle(sprintf('Last Iteration Reconstruction (out of %d iterations)', num_iter), ...
        'FontSize', 13, 'FontWeight', 'bold');

% ── Separate figure for 3rd plot only (num_beats = 11) ────────────────────
b_idx = plot_indices(3);
figure('Name', 'Reconstructed ECG (11 beats)');
hold on;
plot(t, ecg_raw,           'Color',[0.8 0.3 0.3 0.6], 'LineWidth', 0.6);
plot(t, ecg_all(:, b_idx), 'Color',[0.1 0.6 0.2],     'LineWidth', 1.2);
title(sprintf('num\\_beats=%d  |  SNR=%.1f dB', 11, SNR_dB_matrix(end, b_idx)), 'FontSize', 12);
xlabel('Time (s)'); ylabel('Amplitude');
xlim([t(1) t(end)]); grid on; box on;
legend('Noisy ECG', 'Reconstructed ECG', 'Location', 'northeast');

%%  my pan tompkins function based on Rangayyan Book

function [out] = pan_tomp(ECG, Fs, plot_option) 

N = length(ECG);
t = linspace(0, N/Fs, N);

% Low-pass Filter

fn_low = 15;
[B1, A1] = butter(3, fn_low*2/Fs, "low");
ecg_lp = filtfilt(B1, A1, ECG);

%Highpass filter
fn_high = 5;
[B2, A2] = butter(3,fn_high*2/Fs, "high");
ecg_hp = filtfilt(B2, A2, ecg_lp);

%Derivative filter: from (4.14) eq 
B3 = [2 1 0 -1 -2]/8;
ecg_df = filtfilt(B3, 1, ecg_hp);

%Squaring
ecg_sq = ecg_df .^2;

%Integration: from (4.15) eq 
windowSize = 150;   % the book says that 30 is suitable for fs = 200. so I chose 150 for fs = 1000
B4 = (1/windowSize)*ones(1,windowSize);
out = filtfilt(B4, 1, ecg_sq);

%normalize final output
out = out/max(out);

%plot
if plot_option == 1
    figure('Name', 'Pan Tompkin step by step! for ECG')
    subplot(321),       plot(t, ECG),           title('Main ECG')
    subplot(322),       plot(t, ecg_lp),       title('Low Passed')
    subplot(323),       plot(t, ecg_hp),      title('High Passed')
    subplot(324),       plot(t, ecg_df),      title('Deivated')
    subplot(325),       plot(t, ecg_sq),      title('Squared')   , ylabel('time(s)')
    subplot(326),       plot(t, out),            title('Integrated(final result)') , ylabel('time(s)')
end
end
