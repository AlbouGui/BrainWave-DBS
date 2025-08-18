% Author: Albert Guillemette BSc MSc, 02.06.2025
% STFT, MORLET WAVELET and HILBERT ANALYSIS OF SEGMENTED LFP DATA WITH PSD
% AND SPECTROGRAMS.
% This script performs SPECTRAL analysis on X-second epochs, calculates power spectral density,
% creates spectrograms/scalograms.

% Check if segmented_data exists
if ~exist('segmented_data', 'var')
    error('segmented_data not found. Please run the segmentation script first.');
end

%% Parameters
Fs = data.BrainSenseTimeDomain(1).SampleRateInHz; % Sampling rate

% Frequency bands of interest (can be modidifed)
bands.delta = [1, 4];
bands.theta = [4, 8];
bands.alpha = [8, 12];
bands.beta = [12, 30];
bands.gamma = [30, 100];

% STFT parameters
window_length = round(1 * Fs);  % 1-second window
overlap = round(window_length * 0.75);  % 75% overlap
nfft = 2^nextpow2(window_length);  % FFT length (power of 2 for efficiency)

% Morlet wavelet parameters
wavelet_freqs = logspace(log10(1), log10(100), 50);  % 50 frequencies from 1 to 100 Hz
wavelet_cycles = 7;  % Number of cycles for Morlet wavelet (good balance of time-freq resolution)

% Hilbert transform filterbank parameters
filterbank_freqs = logspace(log10(1), log10(100), 50);  % Same frequency grid as wavelet
wavelet_cycles_hilbert = 7;  % To mimic Morlet (constant-Q)

% Hilbert: Low-frequency protections
min_bw_hz = 0.5;   % enforce ≥ this bandwidth at low f (prevents razor-thin bands near DC)
min_fl_hz = 0.5;   % avoid DC proximity for bandpass

% Hilbert: Calculate STFT alignment parameters
stft_win = hamming(window_length, 'periodic');
win_energy = sum(stft_win.^2);                 % window energy U
delta_f = Fs / nfft;                          % frequency resolution
one_sided_factor = 1;                         % set to 2 if your STFT is one-sided & doubled
scale_K = win_energy * nfft * one_sided_factor;

fprintf('===  STFT, SPECTROGRAM, MORLET WAVELET AND HILBERT ANALYSIS ===\n');
fprintf('Sampling rate: %d Hz\n', Fs);
fprintf('STFT window: %.1f seconds (%d samples)\n', window_length/Fs, window_length);
fprintf('Overlap: %.1f%% (%d samples)\n', (overlap/window_length)*100, overlap);
fprintf('FFT length: %d\n', nfft);
fprintf('Frequency resolution: %.2f Hz\n', Fs/nfft);
fprintf('Wavelet frequencies: %.1f - %.1f Hz (%d points)\n', ...
    min(wavelet_freqs), max(wavelet_freqs), length(wavelet_freqs));
fprintf('Wavelet cycles: %d\n', wavelet_cycles);
fprintf(' Hilbert filterbank frequencies: %.1f - %.1f Hz (%d points)\n', ...
    min(filterbank_freqs), max(filterbank_freqs), length(filterbank_freqs));
fprintf('Hilbert bandwidth cycles: %d\n', wavelet_cycles_hilbert);
fprintf(' Hilbert: Using Δf: %.6f Hz | nfft: %d | window energy: %.4f | scale_K: %.4g\n', ...
        delta_f, nfft, win_energy, scale_K);

% Display frequency bands
fprintf('\nFrequency bands:\n');
band_names = fieldnames(bands);
for i = 1:length(band_names)
    fprintf('  %s: %.1f - %.1f Hz\n', band_names{i}, bands.(band_names{i})(1), bands.(band_names{i})(2));
end

%% Initialize results structure
stft_results = struct();
wavelet_results = struct();
hilbert_results = struct();

%% Process each condition (to be modified depending on your # of conditions)
for k = 1:4
    if ~isempty(segmented_data.(['condition_' num2str(k)]))
        fprintf('\n=== Processing Condition %d ===\n', k);
        
        condition_data = segmented_data.(['condition_' num2str(k)]);
        num_epochs = condition_data.num_epochs;
        
        % Initialize arrays for STFT analysis
        all_psds = [];
        all_freqs = [];
        band_powers = struct();
        
        % Initialize arrays for wavelet analysis
        all_scalograms = [];
        wavelet_band_powers = struct();
        
        % Initialize arrays for Hilbert analysis
        epoch_length = size(condition_data.epochs, 2);
        time_vec = (0:epoch_length-1) / Fs;
        all_envelopes = zeros(length(filterbank_freqs), epoch_length, num_epochs);
        all_psds_hilbert = zeros(length(filterbank_freqs), num_epochs);
        hilbert_band_powers = struct();
        
        % ENBW per channel (Hz) for filtfilt chain (|H|^4)
        enbw_vec = zeros(length(filterbank_freqs), 1);
        
        % Initialize band power arrays
        for b = 1:length(band_names)
            band_powers.(band_names{b}) = zeros(num_epochs, 1);
            wavelet_band_powers.(band_names{b}) = zeros(num_epochs, 1);
            hilbert_band_powers.(band_names{b}) = zeros(num_epochs, 1);
        end

        fprintf('Processing %d epochs...\n', num_epochs);
        
        % Process each epoch
        for ep = 1:num_epochs
            epoch_data = condition_data.epochs(ep, :);
            epoch_length = length(epoch_data);
            time_vec = (0:epoch_length-1) / Fs;
            
            %% STFT Analysis
            % Perform STFT
            [S, F, T] = spectrogram(epoch_data, window_length, overlap, nfft, Fs);
            
            % Calculate PSD in V²/Hz (average across time)
            psd = mean(abs(S).^2, 2) * (2 / (Fs * win_energy));  % Normalize by window energy and Fs
            
            % Store STFT results
            if ep == 1
                all_freqs = F;
                all_psds = zeros(length(F), num_epochs);
                % Initialize spectrogram storage
                stft_spectrograms = zeros(length(F), length(T), num_epochs);
                stft_time_vec = T;
            end
            all_psds(:, ep) = psd;
            stft_spectrograms(:, :, ep) = abs(S).^2;
            
                %% Morlet Wavelet Analysis
                % Create Morlet wavelets and compute wavelet transform with ENBW normalization
                wavelet_coefs = zeros(length(wavelet_freqs), epoch_length);
                
                % Initialize ENBW array for the first epoch
                if ep == 1
                    enbw_wav = zeros(length(wavelet_freqs), 1); % Hz
                end
                
                for freq_idx = 1:length(wavelet_freqs)
                    freq = wavelet_freqs(freq_idx);
                    % Create Morlet wavelet
                    sigma = wavelet_cycles / (2 * pi * freq);
                    wavelet_time = -3*sigma:1/Fs:3*sigma;
                    morlet_wavelet = exp(2i * pi * freq * wavelet_time) .* ...
                        exp(-wavelet_time.^2 / (2 * sigma^2));
                    
                    % Normalize wavelet to unit energy (sum(|ψ|^2)*dt = 1)
                    norm_c = sqrt(Fs / sum(abs(morlet_wavelet).^2));
                    morlet_wavelet = morlet_wavelet * norm_c;
                    
                    % Compute ENBW for this wavelet (only for first epoch)
                    if ep == 1
                        enbw_wav(freq_idx) = Fs * sum(abs(morlet_wavelet).^2); % ∫|H(f)|^2 df
                    end
                    
                    % Mirror padding to avoid edge effects
                    pad = min(round(2*Fs), epoch_length-1);
                    xpad = [fliplr(epoch_data(1:pad)), epoch_data, fliplr(epoch_data(end-pad+1:end))];
                    
                    % Convolve on padded signal, then crop to original length
                    conv_pad = conv(xpad, morlet_wavelet, 'same');
                    convolved = conv_pad(pad+1 : pad+epoch_length);
                    wavelet_coefs(freq_idx, :) = convolved;
                end
                
                % Calculate power (magnitude squared)
                wavelet_power = abs(wavelet_coefs).^2; % [F x time]
                
                % Calculate wavelet PSD in V²/Hz (using ENBW)
                wavelet_psd = mean(wavelet_power, 2) ./ enbw_wav; % V²/Hz
                
                % Store wavelet PSD alongside the scalogram data
                if ep == 1
                    all_wavelet_psds = zeros(length(wavelet_freqs), num_epochs);
                end
                all_wavelet_psds(:, ep) = wavelet_psd;
                
                % Store wavelet results
                if ep == 1
                    all_scalograms = zeros(length(wavelet_freqs), epoch_length, num_epochs);
                end
                all_scalograms(:, :, ep) = wavelet_power;
            %% Hilbert Transform Analysis
            % Robust mirror padding: ≥ ~3 cycles of slowest band, ≤ N-1
            pad = round(2*Fs);
            pad = min(pad, numel(epoch_data)-1);
            pad = min(max(pad, round(3*Fs / max(min_fl_hz, min(filterbank_freqs)))), numel(epoch_data)-1);
            sig_padded = [fliplr(epoch_data(1:pad)), epoch_data, fliplr(epoch_data(end-pad+1:end))];
            
            for f_i = 1:length(filterbank_freqs)
                f0 = filterbank_freqs(f_i);

                % Constant-Q geometric edges with LP fallback + min bandwidth clamp
                relBW = 2 / wavelet_cycles_hilbert;                      % ~2/cycles
                r     = 0.5 * (relBW + sqrt(relBW^2 + 4));              % multiplicative half-width
                fl_raw = f0 / r;
                fh_raw = f0 * r;

                % Enforce a minimum absolute bandwidth
                if fh_raw - fl_raw < min_bw_hz
                    mid = (fh_raw + fl_raw)/2;
                    fl_raw = max(mid - min_bw_hz/2, 0);
                    fh_raw = mid + min_bw_hz/2;
                end

                fh = min(fh_raw, Fs/2 - 1);
                if fl_raw < min_fl_hz
                    % Low-pass fallback near DC
                    [b_f, a_f] = butter(4, fh/(Fs/2), 'low');
                    fgrid = linspace(0, fh, 256);
                    lp_mode = true;
                else
                    fl = max(fl_raw, min_fl_hz);
                    if fh <= fl, fh = min(fl + min_bw_hz, Fs/2 - 1); end
                    [b_f, a_f] = butter(4, [fl fh]/(Fs/2), 'bandpass');
                    fgrid = linspace(fl, fh, 256);
                    lp_mode = false;
                end

                % Zero-phase double pass
                filt_padded = filtfilt(b_f, a_f, sig_padded);
                filt_sig    = filt_padded(pad+1 : pad+epoch_length);

                % Analytic envelope POWER
                analytic_sig = hilbert(filt_sig);
                env = abs(analytic_sig).^2;
                all_envelopes(f_i, :, ep) = env;

                % ENBW for this filter (single-pass response integrated; filtfilt => |H|^4)
                if ep == 1
                    Hgrid   = freqz(b_f, a_f, fgrid, Fs);            % single-pass response
                    enbw_hz = trapz(fgrid, abs(Hgrid).^4);           % Hz
                    enbw_vec(f_i) = max(enbw_hz, eps);
                end
            end

            % Calculate Hilbert PSD in V²/Hz
            hilbert_psd = mean(all_envelopes(:,:,ep), 2) ./ enbw_vec;  % V²/Hz
            all_psds_hilbert(:, ep) = hilbert_psd;

            %% Calculate band powers for all 3 methods

            % STFT band power
            for b = 1:length(band_names)
                band_name = band_names{b};
                freq_range = bands.(band_name);
                freq_idx_stft = (F >= freq_range(1)) & (F <= freq_range(2));
                if sum(freq_idx_stft) > 0
                    % Integrate PSD over the frequency band (trapezoidal rule)
                    band_powers.(band_name)(ep) = trapz(F(freq_idx_stft), psd(freq_idx_stft));
                else
                    band_powers.(band_name)(ep) = 0;
                end
            end
            
            % Wavelet band power
            freq_idx_wav = (wavelet_freqs >= freq_range(1)) & (wavelet_freqs <= freq_range(2));
            if sum(freq_idx_wav) > 0
                wavelet_band_powers.(band_name)(ep) = trapz(wavelet_freqs(freq_idx_wav), wavelet_psd(freq_idx_wav));
            else
                wavelet_band_powers.(band_name)(ep) = 0;
            end
            
            % Hilbert band power
            freq_idx_hil = (filterbank_freqs >= freq_range(1)) & (filterbank_freqs <= freq_range(2));
            if sum(freq_idx_hil) > 0
                hilbert_band_powers.(band_name)(ep) = trapz(filterbank_freqs(freq_idx_hil), hilbert_psd(freq_idx_hil));
            else
                hilbert_band_powers.(band_name)(ep) = 0;
            end
        end
        
        %% Store STFT results for this condition
        stft_results.(['condition_' num2str(k)]).psds = all_psds;
        stft_results.(['condition_' num2str(k)]).frequencies = all_freqs;
        stft_results.(['condition_' num2str(k)]).band_powers = band_powers;
        stft_results.(['condition_' num2str(k)]).num_epochs = num_epochs;
        stft_results.(['condition_' num2str(k)]).stft_params = struct(...
            'window_length', window_length, 'overlap', overlap, 'nfft', nfft, 'Fs', Fs);
        stft_results.(['condition_' num2str(k)]).spectrograms = stft_spectrograms;
        stft_results.(['condition_' num2str(k)]).time_vec = stft_time_vec;
        
        % Calculate average PSD and spectrogram across all trials
        mean_psd = mean(all_psds, 2);
        std_psd = std(all_psds, 0, 2);
        mean_spectrogram = mean(stft_spectrograms, 3);
        
        stft_results.(['condition_' num2str(k)]).mean_psd = mean_psd;
        stft_results.(['condition_' num2str(k)]).std_psd = std_psd;
        stft_results.(['condition_' num2str(k)]).mean_spectrogram = mean_spectrogram;
        
        %% Store Wavelet results for this condition
        wavelet_results.(['condition_' num2str(k)]).scalograms = all_scalograms;
        wavelet_results.(['condition_' num2str(k)]).frequencies = wavelet_freqs;
        wavelet_results.(['condition_' num2str(k)]).band_powers = wavelet_band_powers;
        wavelet_results.(['condition_' num2str(k)]).num_epochs = num_epochs;
        wavelet_results.(['condition_' num2str(k)]).wavelet_params = struct(...
            'frequencies', wavelet_freqs, 'cycles', wavelet_cycles);
        wavelet_results.(['condition_' num2str(k)]).time_vec = time_vec;
        wavelet_results.(['condition_' num2str(k)]).psds = all_wavelet_psds;
        wavelet_results.(['condition_' num2str(k)]).mean_psd = mean(all_wavelet_psds, 2);
        wavelet_results.(['condition_' num2str(k)]).std_psd = std(all_wavelet_psds, 0, 2);

        % Calculate average scalogram across all trials
        mean_scalogram = mean(all_scalograms, 3);
        wavelet_results.(['condition_' num2str(k)]).mean_scalogram = mean_scalogram;
        
        %% Store Hilbert results for this condition
        hilbert_results.(['condition_' num2str(k)]).psds = all_psds_hilbert;
        hilbert_results.(['condition_' num2str(k)]).frequencies = filterbank_freqs;
        hilbert_results.(['condition_' num2str(k)]).band_powers = hilbert_band_powers;
        hilbert_results.(['condition_' num2str(k)]).num_epochs = num_epochs;
        hilbert_results.(['condition_' num2str(k)]).hilbert_params = struct(...
            'frequencies', filterbank_freqs, 'cycles', wavelet_cycles_hilbert, ...
            'min_bw_hz', min_bw_hz, 'min_fl_hz', min_fl_hz);
        hilbert_results.(['condition_' num2str(k)]).time_vec = time_vec;
        hilbert_results.(['condition_' num2str(k)]).envelopes = all_envelopes;
        hilbert_results.(['condition_' num2str(k)]).enbw_hz = enbw_vec;                % ENBW (Hz) per channel
        hilbert_results.(['condition_' num2str(k)]).delta_f = delta_f;
        hilbert_results.(['condition_' num2str(k)]).window_energy = win_energy;
        hilbert_results.(['condition_' num2str(k)]).nfft = nfft;
        hilbert_results.(['condition_' num2str(k)]).scale_K = scale_K;
        
        % Calculate average spectrogram for Hilbert method (STFT-aligned units)
        mean_spectrogram_hilbert = bsxfun(@rdivide, mean(all_envelopes,3), enbw_vec) * (delta_f * scale_K);
        hilbert_results.(['condition_' num2str(k)]).mean_spectrogram = mean_spectrogram_hilbert;
        hilbert_results.(['condition_' num2str(k)]).mean_psd = mean(all_psds_hilbert, 2);
        hilbert_results.(['condition_' num2str(k)]).std_psd = std(all_psds_hilbert, 0, 2);
        

        %% Calculate average band powers for all methods
        mean_band_powers = struct();
        std_band_powers = struct();
        mean_wavelet_band_powers = struct();
        std_wavelet_band_powers = struct();
        mean_hilbert_band_powers = struct();
        std_hilbert_band_powers = struct();
        
        for b = 1:length(band_names)
            band_name = band_names{b};
            % STFT band powers
            mean_band_powers.(band_name) = mean(band_powers.(band_name));
            std_band_powers.(band_name) = std(band_powers.(band_name));
            % Wavelet band powers
            mean_wavelet_band_powers.(band_name) = mean(wavelet_band_powers.(band_name));
            std_wavelet_band_powers.(band_name) = std(wavelet_band_powers.(band_name));
            % Hilbert band powers
            mean_hilbert_band_powers.(band_name) = mean(hilbert_band_powers.(band_name));
            std_hilbert_band_powers.(band_name) = std(hilbert_band_powers.(band_name));
        end
        
        stft_results.(['condition_' num2str(k)]).mean_band_powers = mean_band_powers;
        stft_results.(['condition_' num2str(k)]).std_band_powers = std_band_powers;
        wavelet_results.(['condition_' num2str(k)]).mean_band_powers = mean_wavelet_band_powers;
        wavelet_results.(['condition_' num2str(k)]).std_band_powers = std_wavelet_band_powers;
        hilbert_results.(['condition_' num2str(k)]).mean_band_powers = mean_hilbert_band_powers;
        hilbert_results.(['condition_' num2str(k)]).std_band_powers = std_hilbert_band_powers;
        
        fprintf('Condition %d processing completed.\n', k);
        
    else
        fprintf('Condition %d: No data available, skipping...\n', k);
        stft_results.(['condition_' num2str(k)]) = [];
        wavelet_results.(['condition_' num2str(k)]) = [];
        hilbert_results.(['condition_' num2str(k)]) = [];
    end
end
fprintf('\n=== ANALYSIS COMPLETE ===\n');
fprintf('Results stored in:\n');
fprintf('  - stft_results: STFT analysis with spectrograms\n');
fprintf('  - wavelet_results: Morlet wavelet analysis with scalograms\n');
fprintf('  - hilbert_results: Hilbert filterbank analysis (STFT-aligned power units)\n');

%% VISUALIZATION FOR ALL CONDITIONS

fprintf('\n=== Creating Visualizations for All Conditions ===\n');

% Find active conditions
active_conditions = [];
for k = 1:4
    if ~isempty(stft_results.(['condition_' num2str(k)]))
        active_conditions = [active_conditions, k];
    end
end

if isempty(active_conditions)
    fprintf('No active conditions found!\n');
    return;
end

fprintf('Active conditions: %s\n', mat2str(active_conditions));

% Colors for different conditions
colors = {'b', 'c', 'g', 'y'};
line_styles = {'-', '--', '-.', ':'};

%% Plot 1a: Individual PSDs for each condition (separate subplots) - STFT
figure('Name', 'Individual Trial PSDs - All Conditions - STFT', 'Position', [50, 50, 1400, 1000]);

num_active = length(active_conditions);
subplot_rows = ceil(num_active / 2);
subplot_cols = min(2, num_active);

for i = 1:num_active
    k = active_conditions(i);
    condition_results = stft_results.(['condition_' num2str(k)]);
    
    subplot(subplot_rows, subplot_cols, i);
    
    % Plot all individual PSDs in light color
    plot(condition_results.frequencies, 10*log10(condition_results.psds), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    hold on;
    
    % Plot mean PSD in bold
    plot(condition_results.frequencies, 10*log10(condition_results.mean_psd), colors{mod(i-1,4)+1}, 'LineWidth', 2);
    
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title(['Condition ' num2str(k) ' - All Trial PSDs (N=' num2str(condition_results.num_epochs) ')']);
    %legend('Individual Trials', 'Mean PSD', 'Location', 'best');
    grid off;
    xlim([0, 100]);
end

%% Plot 1b: Individual PSDs for each condition (separate subplots)- Wavelet
figure('Name', 'Individual Trial PSDs - All Conditions - Wavelet', 'Position', [50, 50, 1400, 1000]);

num_active = length(active_conditions);
subplot_rows = ceil(num_active / 2);
subplot_cols = min(2, num_active);

for i = 1:num_active
    k = active_conditions(i);
    condition_results = wavelet_results.(['condition_' num2str(k)]);
    
    subplot(subplot_rows, subplot_cols, i);
    
    % Plot all individual PSDs
    plot(condition_results.frequencies, 10*log10(condition_results.psds), ...
         'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    hold on;
    
    % Plot mean PSD in bold 
    plot(condition_results.frequencies, 10*log10(condition_results.mean_psd), ...
         colors{mod(i-1,4)+1}, 'LineWidth', 2);
    
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title(['Condition ' num2str(k) ' - All Trial PSDs (N=' num2str(condition_results.num_epochs) ')']);
    grid off;
    xlim([0, 100]);
end

%% Plot 1c: Individual PSDs for each condition (separate subplots)- Hilbert
figure('Name', 'Individual Trial PSDs - All Conditions - Hilbert', 'Position', [50, 50, 1400, 1000]);

num_active = length(active_conditions);
subplot_rows = ceil(num_active / 2);
subplot_cols = min(2, num_active);

for i = 1:num_active
    k = active_conditions(i);
    condition_results = hilbert_results.(['condition_' num2str(k)]);
    
    subplot(subplot_rows, subplot_cols, i);
    
    % Plot all individual PSDs in light color
    plot(condition_results.frequencies, 10*log10(condition_results.psds), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    hold on;
    
    % Plot mean PSD in bold
    plot(condition_results.frequencies, 10*log10(condition_results.mean_psd), colors{mod(i-1,4)+1}, 'LineWidth', 2);
    
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title(['Condition ' num2str(k) ' - All Trial PSDs (N=' num2str(condition_results.num_epochs) ')']);
    %legend('Individual Trials', 'Mean PSD', 'Location', 'best');
    grid off;
    xlim([0, 100]);
end

%% Plot: Combined PSDs from all three techniques per condition
figure('Name', 'Combined PSDs - STFT, Wavelet, and Hilbert', 'Position', [50, 50, 1400, 1000]);
num_active = length(active_conditions);
subplot_rows = ceil(num_active / 2);
subplot_cols = min(2, num_active);

% Set font to Helvetica and size to 14 for the entire figure
set(groot, 'defaultAxesFontName', 'Helvetica');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultTextFontName', 'Helvetica');
set(groot, 'defaultTextFontSize', 18);

% Define colors for each technique
technique_colors = {
    [0.2, 0.4, 0.8],   % Blue for STFT
    [0.8, 0.2, 0.2],   % Red for Wavelet
    [0.2, 0.8, 0.2]    % Green for Hilbert
};

% First pass: collect all data to determine global y-axis limits
all_power_values = [];
for i = 1:num_active
    k = active_conditions(i);
    
    % Get results from all three techniques
    stft_condition = stft_results.(['condition_' num2str(k)]);
    wavelet_condition = wavelet_results.(['condition_' num2str(k)]);
    hilbert_condition = hilbert_results.(['condition_' num2str(k)]);
    
    % Collect all power values for y-axis scaling
    all_power_values = [all_power_values; 10*log10(stft_condition.mean_psd)];
   wavelet_mean_psd = wavelet_condition.mean_psd;
    all_power_values = [all_power_values; 10*log10(wavelet_mean_psd)];
    all_power_values = [all_power_values; 10*log10(hilbert_condition.mean_psd)];
end

% Calculate global y-axis limits
y_min = min(all_power_values);
y_max = max(all_power_values);
y_range = y_max - y_min;
y_limits = [y_min - 0.05*y_range, y_max + 0.05*y_range]; % Add 5% padding

% Second pass: create plots with consistent y-axis
for i = 1:num_active
    k = active_conditions(i);
    
    % Get results from all three techniques
    stft_condition = stft_results.(['condition_' num2str(k)]);
    wavelet_condition = wavelet_results.(['condition_' num2str(k)]);
    hilbert_condition = hilbert_results.(['condition_' num2str(k)]);
    
    subplot(subplot_rows, subplot_cols, i);
    
    % Plot STFT mean PSD
    plot(stft_condition.frequencies, 10*log10(stft_condition.mean_psd), ...
         'Color', technique_colors{1}, 'LineWidth', 2, 'DisplayName', 'STFT');
    hold on;
    
    % Plot Wavelet mean PSD
    plot(wavelet_condition.frequencies, 10*log10(wavelet_condition.mean_psd), ...
        'Color', technique_colors{2}, 'LineWidth', 2, 'DisplayName', 'Wavelet');
    
    % Plot Hilbert mean PSD
    plot(hilbert_condition.frequencies, 10*log10(hilbert_condition.mean_psd), ...
         'Color', technique_colors{3}, 'LineWidth', 2, 'DisplayName', 'Hilbert');
    
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title(['Condition ' num2str(k) ' - Combined PSDs (N=' num2str(stft_condition.num_epochs) ')']);
    
        % Only show legend on the first subplot
    if i == 1
        legend('Location', 'best');
    end
    
    grid off;
    xlim([0, 100]);
    ylim(y_limits); % Set consistent y-axis limits
    
    hold off;
end

% Add a main title for the entire figure
sgtitle('Power Spectral Density Comparison');

%% Plot 2a: STFT Spectrograms (Average across trials) - ORIGINAL
figure('Name', 'Average STFT Spectrograms - All Conditions', 'Position', [100, 100, 1400, 1000]);

for i = 1:num_active
    k = active_conditions(i);
    condition_results = stft_results.(['condition_' num2str(k)]);
    
    subplot(subplot_rows, subplot_cols, i);
    
    % Plot average spectrogram
    imagesc(condition_results.time_vec, condition_results.frequencies, ...
        10*log10(condition_results.mean_spectrogram));
    axis xy;
    colorbar;
    colormap('parula');
    
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['Condition ' num2str(k) ' - Average STFT Spectrogram']);
    ylim([0, 100]);
    clim([-40, 40]); % Adjust color limits as needed
end

%% Plot 2b: STFT Spectrograms (Average across trials) - RELATIVE POWER
figure('Name', 'Average STFT Spectrograms - All Conditions (RELATIVE POWER)', 'Position', [120, 120, 1400, 1000]);

for i = 1:num_active
    k = active_conditions(i);
    condition_results = stft_results.(['condition_' num2str(k)]);
    
    subplot(subplot_rows, subplot_cols, i);
    
    % Convert to relative power (divide by maximum)
    spectrogram_data = condition_results.mean_spectrogram;
    max_power = max(spectrogram_data(:));
    relative_spectrogram = spectrogram_data / max_power;
    
    % Plot relative spectrogram
    imagesc(condition_results.time_vec, condition_results.frequencies, relative_spectrogram);
    axis xy;
    colorbar;
    colormap('parula');
    
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['Condition ' num2str(k) ' - Average STFT Spectrogram (Relative Power)']);
    ylim([0, 100]);
    clim([0, 1]); % Relative power from 0 to 1
    
    % Add colorbar label
    c = colorbar;
    c.Label.String = 'Relative Power';
end

%% Plot 3a: Morlet Wavelet Scalograms (Average across trials) - ORIGINAL
figure('Name', 'Average Morlet Wavelet Scalograms - All Conditions', 'Position', [150, 150, 1400, 1000]);

for i = 1:num_active
    k = active_conditions(i);
    wavelet_condition_results = wavelet_results.(['condition_' num2str(k)]);
    
    subplot(subplot_rows, subplot_cols, i);
    
    % Plot average scalogram
    imagesc(wavelet_condition_results.time_vec, wavelet_condition_results.frequencies, ...
        10*log10(wavelet_condition_results.mean_scalogram));
    axis xy;
    colorbar;
    colormap('hot');
    
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['Condition ' num2str(k) ' - Average Morlet Wavelet Scalogram']);
    set(gca, 'YScale', 'log'); % Log scale for better visualization
    ylim([1, 100]);
    clim([-40, 40]); % Adjust color limits as needed
end

%% Plot 3b: Morlet Wavelet Scalograms (Average across trials) - RELATIVE POWER
figure('Name', 'Average Morlet Wavelet Scalograms - All Conditions (RELATIVE POWER)', 'Position', [170, 170, 1400, 1000]);

for i = 1:num_active
    k = active_conditions(i);
    wavelet_condition_results = wavelet_results.(['condition_' num2str(k)]);
    
    subplot(subplot_rows, subplot_cols, i);
    
    % Convert to relative power (divide by maximum)
    scalogram_data = wavelet_condition_results.mean_scalogram;
    max_power = max(scalogram_data(:));
    relative_scalogram = scalogram_data / max_power;
    
    % Plot relative scalogram
    imagesc(wavelet_condition_results.time_vec, wavelet_condition_results.frequencies, relative_scalogram);
    axis xy;
    colorbar;
    colormap('hot');
    
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['Condition ' num2str(k) ' - Average Morlet Wavelet Scalogram (Relative Power)']);
    set(gca, 'YScale', 'log'); % Log scale for better visualization
    ylim([1, 100]);
    clim([0, 1]); % Relative power from 0 to 1
    
    % Add colorbar label
    c = colorbar;
    c.Label.String = 'Relative Power';
end

%% Plot 4a: Hilbert Transform Spectrograms (Average across trials) - ORIGINAL
figure('Name', 'Average Hilbert Transform Spectrograms - All Conditions', 'Position', [190, 190, 1400, 1000]);

for i = 1:num_active
    k = active_conditions(i);
    hilbert_condition_results = hilbert_results.(['condition_' num2str(k)]);
    
    subplot(subplot_rows, subplot_cols, i);
    
    % Plot average spectrogram
    imagesc(hilbert_condition_results.time_vec, hilbert_condition_results.frequencies, ...
        10*log10(hilbert_condition_results.mean_spectrogram));
    axis xy;
    colorbar;
    colormap('parula');
    
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['Condition ' num2str(k) ' - Average Hilbert Transform Spectrogram']);
    ylim([0, 100]);
    clim([0, 40]); % Adjust color limits as needed
end

%% Plot 4b: Hilbert Transform Spectrograms (Average across trials) - RELATIVE POWER
figure('Name', 'Average Hilbert Transform Spectrograms - All Conditions (RELATIVE POWER)', 'Position', [210, 210, 1400, 1000]);

for i = 1:num_active
    k = active_conditions(i);
    hilbert_condition_results = hilbert_results.(['condition_' num2str(k)]);
    
    subplot(subplot_rows, subplot_cols, i);
    
    % Convert to relative power (divide by maximum)
    spectrogram_data = hilbert_condition_results.mean_spectrogram;
    max_power = max(spectrogram_data(:));
    relative_spectrogram = spectrogram_data / max_power;
    
    % Plot relative spectrogram
    imagesc(hilbert_condition_results.time_vec, hilbert_condition_results.frequencies, relative_spectrogram);
    axis xy;
    colorbar;
    colormap('parula');
    
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['Condition ' num2str(k) ' - Average Hilbert Transform Spectrogram (Relative Power)']);
    ylim([0, 100]);
    clim([0, 1]); % Relative power from 0 to 1
    
    % Add colorbar label
    c = colorbar;
    c.Label.String = 'Relative Power';
end

%% Interactive Condition Selection for Individual Trial Spectrograms - ORIGINAL - STFT
% Allow user to choose which condition to display for individual trials
fprintf('\n=== Individual Trial Spectrograms (ORIGINAL) ===\n');
fprintf('Available conditions: %s\n', mat2str(active_conditions));

% User selection for condition (with default to first active condition)
selected_condition = active_conditions(1); % Default
user_input = input(sprintf('Enter condition number to display individual trial spectrograms (STFT) (default: %d): ', selected_condition), 's');
if ~isempty(user_input)
    temp_condition = str2double(user_input);
    if ismember(temp_condition, active_conditions)
        selected_condition = temp_condition;
    else
        fprintf('Invalid condition selected. Using default condition %d.\n', selected_condition);
    end
end

% Get the results for selected condition
condition_results = stft_results.(['condition_' num2str(selected_condition)]);

% Display ALL trials
% Calculate subplot arrangement for all trials
num_trials = condition_results.num_epochs;
subplot_cols_trials = ceil(sqrt(num_trials));
subplot_rows_trials = ceil(num_trials / subplot_cols_trials);

% Create figure for all individual trial spectrograms
figure('Name', ['All Individual Trial Spectrograms - Condition ' num2str(selected_condition)], ...
    'Position', [250, 250, 1600, 1200]);

for trial = 1:num_trials
    subplot(subplot_rows_trials, subplot_cols_trials, trial);
    
    imagesc(condition_results.time_vec, condition_results.frequencies, ...
        10*log10(condition_results.spectrograms(:, :, trial)));
    axis xy;
    colorbar;
    colormap('parula');
    
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['Trial ' num2str(trial)]);
    ylim([0, 100]);
    clim([-40, 40]);
end

% Add a main title for the entire figure
sgtitle(['Condition ' num2str(selected_condition) ' - All Individual Trial STFT Spectrograms (N=' num2str(num_trials) ')']);

%% Individual Trial Spectrograms - RELATIVE POWER
fprintf('\n=== Individual Trial Spectrograms (RELATIVE POWER) ===\n');

% Calculate global maximum across ALL trials for consistent normalization
global_max_stft = max(condition_results.spectrograms(:));
fprintf('Global STFT maximum power: %.2e\n', global_max_stft);

% Create figure for all individual trial spectrograms with relative power
figure('Name', ['All Individual Trial Spectrograms - Condition ' num2str(selected_condition) ' (RELATIVE POWER)'], ...
    'Position', [270, 270, 1600, 1200]);

for trial = 1:num_trials
    subplot(subplot_rows_trials, subplot_cols_trials, trial);
    
    % Convert to relative power using global maximum
    trial_spectrogram = condition_results.spectrograms(:, :, trial);
    relative_trial_spectrogram = trial_spectrogram / global_max_stft;
    
    imagesc(condition_results.time_vec, condition_results.frequencies, relative_trial_spectrogram);
    axis xy;
    colorbar;
    colormap('parula');
    
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['Trial ' num2str(trial)]);
    ylim([0, 100]);
    clim([0, 1]); % Consistent relative power scale
    
    % Add colorbar label for first subplot
    if trial == 1
        c = colorbar;
        c.Label.String = 'Relative Power';
    end
end

% Add a main title for the entire figure
sgtitle(['Condition ' num2str(selected_condition) ' - All Individual Trial STFT Spectrograms - RELATIVE POWER (N=' num2str(num_trials) ')']);

%% Individual Trial Scalograms (wavelet) - ORIGINAL
fprintf('\n=== Individual Trial Scalograms (ORIGINAL) ===\n');
fprintf('Available conditions: %s\n', mat2str(active_conditions));

% User selection for condition (with default to first active condition)
selected_condition_wav = active_conditions(1); % Default
user_input_wav = input(sprintf('Enter condition number to display individual trial scalograms - Wavelet (default: %d): ', selected_condition_wav), 's');
if ~isempty(user_input_wav)
    temp_condition_wav = str2double(user_input_wav);
    if ismember(temp_condition_wav, active_conditions)
        selected_condition_wav = temp_condition_wav;
    else
        fprintf('Invalid condition selected. Using default condition %d.\n', selected_condition_wav);
    end
end

% Get the wavelet results for selected condition
wavelet_condition_results = wavelet_results.(['condition_' num2str(selected_condition_wav)]);

% Calculate consistent color limits across all trials
all_scalograms_db = 10*log10(wavelet_condition_results.scalograms);
global_min = prctile(all_scalograms_db(:), 5);  % 5th percentile to avoid extreme outliers
global_max = prctile(all_scalograms_db(:), 95); % 95th percentile

fprintf('Global color range: %.1f to %.1f dB\n', global_min, global_max);

% Display ALL trials for wavelet scalograms
num_trials_wav = wavelet_condition_results.num_epochs;
subplot_cols_trials_wav = ceil(sqrt(num_trials_wav));
subplot_rows_trials_wav = ceil(num_trials_wav / subplot_cols_trials_wav);

% Create figure for all individual trial scalograms
figure('Name', ['All Individual Trial Scalograms - Condition ' num2str(selected_condition_wav)], ...
    'Position', [300, 300, 1600, 1200]);

for trial = 1:num_trials_wav
    subplot(subplot_rows_trials_wav, subplot_cols_trials_wav, trial);
    
    % Convert to dB for visualization
    scalogram_db = 10*log10(wavelet_condition_results.scalograms(:, :, trial));
    
    imagesc(wavelet_condition_results.time_vec, wavelet_condition_results.frequencies, scalogram_db);
    axis xy;
    colorbar;
    colormap('hot'); % Dark red (low) to white (high)
    
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['Trial ' num2str(trial)]);
    set(gca, 'YScale', 'log');
    ylim([1, 100]);
    
    % Use consistent color limits across all trials
    clim([global_min, global_max]);
end

% Add a main title for the entire figure
sgtitle(['Condition ' num2str(selected_condition_wav) ' - All Individual Trial Morlet Wavelet Scalograms (N=' num2str(num_trials_wav) ')']);

%% Individual Trial Scalograms (wavelet) - RELATIVE POWER
fprintf('\n=== Individual Trial Scalograms (RELATIVE POWER) ===\n');

% Calculate global maximum across ALL scalogram trials for consistent normalization
global_max_wavelet = max(wavelet_condition_results.scalograms(:));
fprintf('Global Wavelet maximum power: %.2e\n', global_max_wavelet);

% Create figure for all individual trial scalograms with relative power
figure('Name', ['All Individual Trial Scalograms - Condition ' num2str(selected_condition_wav) ' (RELATIVE POWER)'], ...
    'Position', [320, 320, 1600, 1200]);

for trial = 1:num_trials_wav
    subplot(subplot_rows_trials_wav, subplot_cols_trials_wav, trial);
    
    % Convert to relative power using global maximum
    trial_scalogram = wavelet_condition_results.scalograms(:, :, trial);
    relative_trial_scalogram = trial_scalogram / global_max_wavelet;
    
    imagesc(wavelet_condition_results.time_vec, wavelet_condition_results.frequencies, relative_trial_scalogram);
    axis xy;
    colorbar;
    colormap('hot'); % Dark red (low) to white (high)
    
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['Trial ' num2str(trial)]);
    set(gca, 'YScale', 'log');
    ylim([1, 100]);
    
    % Use consistent relative power limits across all trials
    clim([0, 1]);
    
    % Add colorbar label for first subplot
    if trial == 1
        c = colorbar;
        c.Label.String = 'Relative Power';
    end
end

% Add a main title for the entire figure
sgtitle(['Condition ' num2str(selected_condition_wav) ' - All Individual Trial Morlet Wavelet Scalograms - RELATIVE POWER (N=' num2str(num_trials_wav) ')']);

%% Individual Trial Spectrogram (Hilbert) - ORIGINAL
fprintf('\n=== Individual Trial Spectrogram (ORIGINAL) ===\n');
fprintf('Available conditions: %s\n', mat2str(active_conditions));

% User selection for condition (with default to first active condition)
selected_condition_wav = active_conditions(1); % Default
user_input_wav = input(sprintf('Enter condition number to display individual trial spectrogram - Hilbert (default: %d): ', selected_condition_wav), 's');
if ~isempty(user_input_wav)
    temp_condition_wav = str2double(user_input_wav);
    if ismember(temp_condition_wav, active_conditions)
        selected_condition_wav = temp_condition_wav;
    else
        fprintf('Invalid condition selected. Using default condition %d.\n', selected_condition_wav);
    end
end

% Get the results for selected condition
condition_results = hilbert_results.(['condition_' num2str(selected_condition)]);

% Display ALL trials
% Calculate subplot arrangement for all trials
num_trials = condition_results.num_epochs;
subplot_cols_trials = ceil(sqrt(num_trials));
subplot_rows_trials = ceil(num_trials / subplot_cols_trials);

% Create figure for all individual trial spectrograms
figure('Name', ['All Individual Trial Spectrograms - Condition ' num2str(selected_condition)], ...
    'Position', [250, 250, 1600, 1200]);


for trial = 1:num_trials
    subplot(subplot_rows_trials, subplot_cols_trials, trial);
    
    trial_spectrogram = condition_results.envelopes(:, :, trial) ./ condition_results.enbw_hz * (condition_results.delta_f * condition_results.scale_K);
    imagesc(condition_results.time_vec, condition_results.frequencies, ...
        10*log10(trial_spectrogram));

clim([-40, 40]);

    axis xy;
    colorbar;
    colormap('parula');
    
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['Trial ' num2str(trial)]);
    ylim([0, 100]);
    clim([0, 40]);
end

%% Individual Trial Spectrograms (Hilbert) - RELATIVE POWER
fprintf('\n=== Individual Trial Spectrograms/Hilbert (RELATIVE POWER) ===\n');

% Calculate STFT-aligned spectrograms for all trials to find global maximum
all_trials_spectrograms = zeros(size(condition_results.envelopes));
for trial = 1:num_trials
    all_trials_spectrograms(:, :, trial) = condition_results.envelopes(:, :, trial) ./ condition_results.enbw_hz * (condition_results.delta_f * condition_results.scale_K);
end

% Calculate global maximum across all trials in STFT-aligned units
global_max_hilbert = max(all_trials_spectrograms(isfinite(all_trials_spectrograms(:))));
if isempty(global_max_hilbert) || ~isfinite(global_max_hilbert)
    warning('Invalid values in Hilbert spectrograms for condition %d. Using default max value.', selected_condition);
    global_max_hilbert = 1; % Fallback value
end
fprintf('Global Hilbert maximum power (STFT-aligned units): %.2e\n', global_max_hilbert);

% Create figure for all individual trial spectrograms with relative power
figure('Name', ['All Individual Trial Spectrograms - Condition ' num2str(selected_condition) ' (RELATIVE POWER - Hilbert)'], ...
    'Position', [270, 270, 1600, 1200]);

for trial = 1:num_trials
    subplot(subplot_rows_trials, subplot_cols_trials, trial);
    
    % Use precomputed STFT-aligned spectrogram for this trial
    trial_spectrogram = all_trials_spectrograms(:, :, trial);
    relative_trial_spectrogram = trial_spectrogram / global_max_hilbert;
    
    imagesc(condition_results.time_vec, condition_results.frequencies, relative_trial_spectrogram);
    axis xy;
    colorbar;
    colormap('parula');
    
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['Trial ' num2str(trial)]);
    ylim([0, 100]);
    clim([0, 1]); % Consistent relative power scale
    
    % Add colorbar label for first subplot
    if trial == 1
        c = colorbar;
        c.Label.String = 'Relative Power';
    end
end

% Add a main title for the entire figure
sgtitle(['Condition ' num2str(selected_condition) ' - All Individual Trial Hilbert Spectrograms - RELATIVE POWER (N=' num2str(num_trials) ')']);%% End
fprintf('\nBoth original and relative power visualizations completed!\n');
fprintf('\nAll visualizations completed with BOTH original and relative power versions!\n');
fprintf('Created enhanced figures with:\n');
fprintf('- Original absolute power visualizations (dB scale)\n');
fprintf('- Relative power normalizations (0-1 scale) for optimal visualization\n');
fprintf('- Consistent color scales within each visualization type\n');
fprintf('- Interactive condition selection for individual trials\n');
fprintf('- Display of ALL trials \n');
fprintf('- Properly arranged subplots for any number of trials\n');
fprintf('- Global normalization ensures fair comparison across trials in relative power plots\n');
fprintf('Results stored in stft_results and wavelet_results structures\n');

% :)