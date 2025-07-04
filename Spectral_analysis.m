% Author: Albert Guillemette BSc MSc, 02.06.2025
% STFT, PSD, SPECTROGRAM AND MORLET WAVELET ANALYSIS OF SEGMENTED LFP DATA
% This script performs STFT on 5-second epochs, calculates power spectral density,
% creates spectrograms, and performs Morlet wavelet analysis with scalograms
% ENHANCED VERSION: Uses relative power (normalized by max) for optimal visualization

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
window_length = round(2 * Fs);  % 2-second window
overlap = round(window_length * 0.75);  % 75% overlap
nfft = 2^nextpow2(window_length);  % FFT length (power of 2 for efficiency)

% Morlet wavelet parameters
wavelet_freqs = logspace(log10(1), log10(100), 50);  % 50 frequencies from 1 to 100 Hz
wavelet_cycles = 7;  % Number of cycles for Morlet wavelet (good balance of time-freq resolution)

fprintf('=== ENHANCED STFT, SPECTROGRAM AND MORLET WAVELET ANALYSIS ===\n');
fprintf('Sampling rate: %d Hz\n', Fs);
fprintf('STFT window: %.1f seconds (%d samples)\n', window_length/Fs, window_length);
fprintf('Overlap: %.1f%% (%d samples)\n', (overlap/window_length)*100, overlap);
fprintf('FFT length: %d\n', nfft);
fprintf('Frequency resolution: %.2f Hz\n', Fs/nfft);
fprintf('Wavelet frequencies: %.1f - %.1f Hz (%d points)\n', ...
    min(wavelet_freqs), max(wavelet_freqs), length(wavelet_freqs));
fprintf('Wavelet cycles: %d\n', wavelet_cycles);

% Display frequency bands
fprintf('\nFrequency bands:\n');
band_names = fieldnames(bands);
for i = 1:length(band_names)
    fprintf('  %s: %.1f - %.1f Hz\n', band_names{i}, bands.(band_names{i})(1), bands.(band_names{i})(2));
end

%% Initialize results structure
stft_results = struct();
wavelet_results = struct();

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
        
        % Initialize band power arrays
        for b = 1:length(band_names)
            band_powers.(band_names{b}) = zeros(num_epochs, 1);
            wavelet_band_powers.(band_names{b}) = zeros(num_epochs, 1);
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
            
            % Calculate Power Spectral Density (average across time)
            psd = mean(abs(S).^2, 2);  % Average power across time windows
            
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
            % Create Morlet wavelets and compute wavelet transform
            wavelet_coefs = zeros(length(wavelet_freqs), epoch_length);
            
            for freq_idx = 1:length(wavelet_freqs)
                freq = wavelet_freqs(freq_idx);
                % Create Morlet wavelet
                sigma = wavelet_cycles / (2 * pi * freq);
                wavelet_time = -3*sigma:1/Fs:3*sigma;
                morlet_wavelet = exp(2i * pi * freq * wavelet_time) .* ...
                    exp(-wavelet_time.^2 / (2 * sigma^2));
                
                % Convolve with signal (using convolution)
                convolved = conv(epoch_data, morlet_wavelet, 'same');
                wavelet_coefs(freq_idx, :) = convolved;
            end
            
            % Calculate power (magnitude squared)
            wavelet_power = abs(wavelet_coefs).^2;
            
            % Store wavelet results
            if ep == 1
                all_scalograms = zeros(length(wavelet_freqs), epoch_length, num_epochs);
            end
            all_scalograms(:, :, ep) = wavelet_power;
            
            %% Calculate band powers for both methods
            % STFT band powers
            for b = 1:length(band_names)
                band_name = band_names{b};
                freq_range = bands.(band_name);
                
                % STFT band power
                freq_idx_stft = (F >= freq_range(1)) & (F <= freq_range(2));
                if sum(freq_idx_stft) > 0
                    band_powers.(band_name)(ep) = mean(psd(freq_idx_stft));
                else
                    band_powers.(band_name)(ep) = 0;
                end
                
                % Wavelet band power
                freq_idx_wav = (wavelet_freqs >= freq_range(1)) & (wavelet_freqs <= freq_range(2));
                if sum(freq_idx_wav) > 0
                    wavelet_band_powers.(band_name)(ep) = mean(mean(wavelet_power(freq_idx_wav, :)));
                else
                    wavelet_band_powers.(band_name)(ep) = 0;
                end
            end
            
            if mod(ep, 10) == 0 || ep == num_epochs
                fprintf('  Processed %d/%d epochs\n', ep, num_epochs);
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
        
        % Calculate average scalogram across all trials
        mean_scalogram = mean(all_scalograms, 3);
        wavelet_results.(['condition_' num2str(k)]).mean_scalogram = mean_scalogram;
        
        % Calculate average band powers for both methods
        mean_band_powers = struct();
        std_band_powers = struct();
        mean_wavelet_band_powers = struct();
        std_wavelet_band_powers = struct();
        
        for b = 1:length(band_names)
            band_name = band_names{b};
            % STFT band powers
            mean_band_powers.(band_name) = mean(band_powers.(band_name));
            std_band_powers.(band_name) = std(band_powers.(band_name));
            % Wavelet band powers
            mean_wavelet_band_powers.(band_name) = mean(wavelet_band_powers.(band_name));
            std_wavelet_band_powers.(band_name) = std(wavelet_band_powers.(band_name));
        end
        
        stft_results.(['condition_' num2str(k)]).mean_band_powers = mean_band_powers;
        stft_results.(['condition_' num2str(k)]).std_band_powers = std_band_powers;
        wavelet_results.(['condition_' num2str(k)]).mean_band_powers = mean_wavelet_band_powers;
        wavelet_results.(['condition_' num2str(k)]).std_band_powers = std_wavelet_band_powers;
        
        fprintf('Condition %d processing completed.\n', k);
        
    else
        fprintf('Condition %d: No data available, skipping...\n', k);
        stft_results.(['condition_' num2str(k)]) = [];
        wavelet_results.(['condition_' num2str(k)]) = [];
    end
end

%% ENHANCED VISUALIZATION FOR ALL CONDITIONS

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
colors = {'b', 'r', 'g', 'm'};
line_styles = {'-', '--', '-.', ':'};

%% Plot 1: Individual PSDs for each condition (separate subplots)
figure('Name', 'Individual Trial PSDs - All Conditions', 'Position', [50, 50, 1400, 1000]);

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

%% Plot 2: STFT Spectrograms (Average across trials) - ORIGINAL
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

%% Plot 3: Morlet Wavelet Scalograms (Average across trials) - ORIGINAL
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
    clim([-60, 0]); % Adjust color limits as needed
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

%% Interactive Condition Selection for Individual Trial Spectrograms - ORIGINAL
% Allow user to choose which condition to display for individual trials
fprintf('\n=== Individual Trial Spectrograms (ORIGINAL) ===\n');
fprintf('Available conditions: %s\n', mat2str(active_conditions));

% User selection for condition (with default to first active condition)
selected_condition = active_conditions(1); % Default
user_input = input(sprintf('Enter condition number to display individual trial spectrograms (default: %d): ', selected_condition), 's');
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

%% Individual Trial Scalograms - ORIGINAL
fprintf('\n=== Individual Trial Scalograms (ORIGINAL) ===\n');
fprintf('Available conditions: %s\n', mat2str(active_conditions));

% User selection for condition (with default to first active condition)
selected_condition_wav = active_conditions(1); % Default
user_input_wav = input(sprintf('Enter condition number to display individual trial scalograms (default: %d): ', selected_condition_wav), 's');
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

%% Individual Trial Scalograms - RELATIVE POWER
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
