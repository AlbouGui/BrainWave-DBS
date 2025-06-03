%% STFT AND PSD ANALYSIS OF SEGMENTED LFP DATA
% This script performs STFT on 5-second epochs and calculates power spectral density
% Calculates average power in predefined frequency bands

% Check if segmented_data exists
if ~exist('segmented_data', 'var')
    error('segmented_data not found. Please run the segmentation script first.');
end

%% Parameters
Fs = data.BrainSenseTimeDomain(1).SampleRateInHz; % Sampling rate

% Frequency bands of interest
bands.delta = [1, 4];
bands.theta = [4, 8];
bands.alpha = [8, 12];
bands.beta = [12, 30];
bands.gamma = [30, 100];

% STFT parameters
window_length = round(2 * Fs);  % 2-second window
overlap = round(window_length * 0.75);  % 75% overlap
nfft = 2^nextpow2(window_length);  % FFT length (power of 2 for efficiency)

fprintf('=== STFT AND PSD ANALYSIS ===\n');
fprintf('Sampling rate: %d Hz\n', Fs);
fprintf('STFT window: %.1f seconds (%d samples)\n', window_length/Fs, window_length);
fprintf('Overlap: %.1f%% (%d samples)\n', (overlap/window_length)*100, overlap);
fprintf('FFT length: %d\n', nfft);
fprintf('Frequency resolution: %.2f Hz\n', Fs/nfft);

% Display frequency bands
fprintf('\nFrequency bands:\n');
band_names = fieldnames(bands);
for i = 1:length(band_names)
    fprintf('  %s: %.1f - %.1f Hz\n', band_names{i}, bands.(band_names{i})(1), bands.(band_names{i})(2));
end

%% Initialize results structure
stft_results = struct();

%% Process each channel
for k = 1:4
    if ~isempty(segmented_data.(['channel_' num2str(k)]))
        fprintf('\n=== Processing Channel %d ===\n', k);
        
        channel_data = segmented_data.(['channel_' num2str(k)]);
        num_epochs = channel_data.num_epochs;
        
        % Initialize arrays for this channel
        all_psds = [];
        all_freqs = [];
        band_powers = struct();
        
        % Initialize band power arrays
        for b = 1:length(band_names)
            band_powers.(band_names{b}) = zeros(num_epochs, 1);
        end
        
        fprintf('Processing %d epochs...\n', num_epochs);
        
        % Process each epoch
        for ep = 1:num_epochs
            epoch_data = channel_data.epochs(ep, :);
            
            % Perform STFT
            [S, F, T] = spectrogram(epoch_data, window_length, overlap, nfft, Fs);
            
            % Calculate Power Spectral Density (average across time)
            psd = mean(abs(S).^2, 2);  % Average power across time windows
            
            % Store results
            if ep == 1
                all_freqs = F;
                all_psds = zeros(length(F), num_epochs);
            end
            all_psds(:, ep) = psd;
            
            % Calculate band powers for this epoch
            for b = 1:length(band_names)
                band_name = band_names{b};
                freq_range = bands.(band_name);
                
                % Find frequency indices within the band
                freq_idx = (F >= freq_range(1)) & (F <= freq_range(2));
                
                % Calculate average power in this band
                if sum(freq_idx) > 0
                    band_powers.(band_name)(ep) = mean(psd(freq_idx));
                else
                    band_powers.(band_name)(ep) = 0;
                    if ep == 1
                        fprintf('Warning: No frequencies found in %s band (%.1f-%.1f Hz)\n', ...
                            band_name, freq_range(1), freq_range(2));
                    end
                end
            end
            
            if mod(ep, 10) == 0 || ep == num_epochs
                fprintf('  Processed %d/%d epochs\n', ep, num_epochs);
            end
        end
        
        % Store results for this channel
        stft_results.(['channel_' num2str(k)]).psds = all_psds;
        stft_results.(['channel_' num2str(k)]).frequencies = all_freqs;
        stft_results.(['channel_' num2str(k)]).band_powers = band_powers;
        stft_results.(['channel_' num2str(k)]).num_epochs = num_epochs;
        stft_results.(['channel_' num2str(k)]).stft_params = struct(...
            'window_length', window_length, 'overlap', overlap, 'nfft', nfft, 'Fs', Fs);
        
        % Calculate average PSD across all trials
        mean_psd = mean(all_psds, 2);
        std_psd = std(all_psds, 0, 2);
        stft_results.(['channel_' num2str(k)]).mean_psd = mean_psd;
        stft_results.(['channel_' num2str(k)]).std_psd = std_psd;
        
        % Calculate average band powers
        mean_band_powers = struct();
        std_band_powers = struct();
        for b = 1:length(band_names)
            band_name = band_names{b};
            mean_band_powers.(band_name) = mean(band_powers.(band_name));
            std_band_powers.(band_name) = std(band_powers.(band_name));
        end
        stft_results.(['channel_' num2str(k)]).mean_band_powers = mean_band_powers;
        stft_results.(['channel_' num2str(k)]).std_band_powers = std_band_powers;
        
        fprintf('Channel %d processing completed.\n', k);
        
    else
        fprintf('Channel %d: No data available, skipping...\n', k);
        stft_results.(['channel_' num2str(k)]) = [];
    end
end

%% ENHANCED VISUALIZATION FOR ALL CHANNELS

fprintf('\n=== Creating Visualizations for All Channels ===\n');

% Find active channels
active_channels = [];
for k = 1:4
    if ~isempty(stft_results.(['channel_' num2str(k)]))
        active_channels = [active_channels, k];
    end
end

if isempty(active_channels)
    fprintf('No active channels found!\n');
    return;
end

fprintf('Active channels: %s\n', mat2str(active_channels));

% Colors for different channels
colors = {'b', 'r', 'g', 'm'};
line_styles = {'-', '--', '-.', ':'};

%% Plot 1: Individual PSDs for each channel (separate subplots)
figure('Name', 'Individual Trial PSDs - All Channels', 'Position', [50, 50, 1400, 1000]);

num_active = length(active_channels);
subplot_rows = ceil(num_active / 2);
subplot_cols = min(2, num_active);

for i = 1:num_active
    k = active_channels(i);
    channel_results = stft_results.(['channel_' num2str(k)]);
    
    subplot(subplot_rows, subplot_cols, i);
    
    % Plot all individual PSDs in light color
    plot(channel_results.frequencies, 10*log10(channel_results.psds), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    hold on;
    
    % Plot mean PSD in bold
    plot(channel_results.frequencies, 10*log10(channel_results.mean_psd), colors{mod(i-1,4)+1}, 'LineWidth', 2);
    
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title(['Channel ' num2str(k) ' - All Trial PSDs (N=' num2str(channel_results.num_epochs) ')']);
    legend('Individual Trials', 'Mean PSD', 'Location', 'best');
    grid on;
    xlim([0, 100]);
end

%% Plot 2: Mean PSDs with error bars for each channel
figure('Name', 'Mean PSDs with Error Bars - All Channels', 'Position', [100, 100, 1400, 1000]);

for i = 1:num_active
    k = active_channels(i);
    channel_results = stft_results.(['channel_' num2str(k)]);
    
    subplot(subplot_rows, subplot_cols, i);
    
    freq_vec = channel_results.frequencies;
    mean_psd_db = 10*log10(channel_results.mean_psd);
    std_psd_db = 10*log10(channel_results.std_psd);
    
    % Create error band
    color_light = [0.8 0.8 1];
    if i == 2, color_light = [1 0.8 0.8]; end
    if i == 3, color_light = [0.8 1 0.8]; end
    if i == 4, color_light = [1 0.8 1]; end
    
    fill([freq_vec; flipud(freq_vec)], ...
         [mean_psd_db + std_psd_db; flipud(mean_psd_db - std_psd_db)], ...
         color_light, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    hold on;
    
    % Plot mean line
    plot(freq_vec, mean_psd_db, colors{mod(i-1,4)+1}, 'LineWidth', 2);
    
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title(['Channel ' num2str(k) ' - Mean PSD ± STD']);
    grid on;
    xlim([0, 100]);
end

%% Plot 3: Band powers across trials for each channel
figure('Name', 'Band Powers Across Trials - All Channels', 'Position', [150, 150, 1400, 1000]);

for i = 1:num_active
    k = active_channels(i);
    channel_results = stft_results.(['channel_' num2str(k)]);
    
    subplot(subplot_rows, subplot_cols, i);
    
    % Prepare band power data for boxplot
    band_data = [];
    band_labels = {};
    for b = 1:length(band_names)
        band_data = [band_data, channel_results.band_powers.(band_names{b})];
        band_labels{b} = band_names{b};
    end
    
    boxplot(band_data, 'Labels', band_labels);
    ylabel('Power');
    title(['Channel ' num2str(k) ' - Band Powers Across Trials']);
    xtickangle(45);
    grid on;
end

%% Plot 4: Mean band powers comparison for each channel
figure('Name', 'Mean Band Powers - All Channels', 'Position', [200, 200, 1400, 1000]);

for i = 1:num_active
    k = active_channels(i);
    channel_results = stft_results.(['channel_' num2str(k)]);
    
    subplot(subplot_rows, subplot_cols, i);
    
    % Get mean and std for each band
    mean_powers = [];
    std_powers = [];
    for b = 1:length(band_names)
        mean_powers(b) = channel_results.mean_band_powers.(band_names{b});
        std_powers(b) = channel_results.std_band_powers.(band_names{b});
    end
    
    % Create bar plot with error bars
    bar(1:length(band_names), mean_powers, 'FaceColor', colors{mod(i-1,4)+1}, 'FaceAlpha', 0.7);
    hold on;
    errorbar(1:length(band_names), mean_powers, std_powers, 'k.', 'LineWidth', 1.5);
    
    set(gca, 'XTickLabel', band_names);
    ylabel('Average Power');
    title(['Channel ' num2str(k) ' - Mean Band Powers']);
    xtickangle(45);
    grid on;
end

fprintf('All visualizations completed!\n');
fprintf('Created %d figures for %d active channels\n', ...
        3 + (length(active_channels) > 1), length(active_channels));