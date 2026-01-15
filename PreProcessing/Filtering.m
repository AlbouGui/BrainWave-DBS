% Author: Albert Guillemette BSc MSc, 02.06.2025
% FILTERING apply low and high pass butterworth filters, notch filter and
% cardiac artefact management (if required)

% Get sampling rate
Fs = data.BrainSenseTimeDomain(1).SampleRateInHz;
fprintf('Sampling rate: %d Hz\n', Fs);

% Initialize filtered data structures
d.no_60Hz = cell(4,1);
d.butter_filtered = cell(4,1);
d.cleaned_LFP = cell(4,1);  % Added initialization
d.smoothed_final = cell(4,1);  % Renamed for clarity

%% 1) Notch filter at 60 Hz (Electrical interference in North America)
design_filter = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',Fs);

fprintf('Applying 60Hz notch filter...\n');
for k = 1:4
    if d.raw{k}.initialized && ~isempty(d.raw{k}.LFP)
        % Apply notch filter to LFP signal only
        filtered_lfp = filtfilt(design_filter, d.raw{k}.LFP);
        % Store as [LFP; triggers] format for compatibility
        d.no_60Hz{k} = [filtered_lfp; d.raw{k}.triggers];
    else
        d.no_60Hz{k} = [];
    end
end

%% 2) Bandpass filter
% Designing a Butterworth filter (4th order, bandpass 1-100 Hz)
butterworth_filter = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', 1, 'HalfPowerFrequency2', 100, ...
    'DesignMethod', 'butter', 'SampleRate', Fs);

fprintf('Applying bandpass filter (1-100 Hz)...\n');
for k = 1:4
    if ~isempty(d.no_60Hz{k})
        % Apply bandpass filter
        filtered_lfp = filtfilt(butterworth_filter, d.no_60Hz{k}(1,:));
        % Keep original triggers
        d.butter_filtered{k} = [filtered_lfp; d.raw{k}.triggers];
    else
        d.butter_filtered{k} = [];
    end
end

%% 3) Cardiac artifact removal (selective by condition i.e., LeftSTNOff ...)
% Specify which condition need cardiac cleaning based on visual inspection
% Set to empty array [] if no cardiac cleaning needed
condition_needing_cleaning = [4];  % e.g., [1, 3] for conditions 1 and 3
% condition_needing_cleaning = [1, 3];  % Example: clean conditions 1 and 3

% Initialize cleaned_LFP with butter_filtered data
d.cleaned_LFP = d.butter_filtered;

if ~isempty(condition_needing_cleaning)
    fprintf('Applying cardiac artifact cleaning to condition: %s\n', num2str(condition_needing_cleaning));
    
    for k = condition_needing_cleaning  % Only process specified condition
        if ~isempty(d.butter_filtered{k})
            LFP_signal = d.butter_filtered{k}(1,:); % Extract LFP
            
            % Modified Z-score normalization of the LFP signal
            median_signal = median(LFP_signal);
            mad_signal = mad(abs(LFP_signal - median_signal));
            
            % Avoid division by zero
            if mad_signal == 0
                mad_signal = std(LFP_signal);
            end
            
            modified_z_scored_signal = (LFP_signal - median_signal) / mad_signal;
            
            % Parameters for peak detection 
            min_peak_height = 2.5; % Now in terms of MAD units, can be modified upon peaks visualisation
            min_peak_distance = round(0.5 * Fs); % 0.5 seconds between peaks; to be modified depending on the participant HR
            
            % Detect positive and negative peaks
            [peaks_positive, locs_positive] = findpeaks(modified_z_scored_signal, ...
                'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_distance);
            [peaks_negative, locs_negative] = findpeaks(-modified_z_scored_signal, ...
                'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_distance);
            
            % Determine orientation based on peak strength
            if length(locs_positive) >= length(locs_negative)
                detected_peaks = locs_positive;
                orientation = 'positive';
            else
                detected_peaks = locs_negative;
                orientation = 'negative';
            end
            
            % Only proceed if we have enough peaks
            if length(detected_peaks) >= 3
                % Plot detected peaks
                figure;
                plot(modified_z_scored_signal);
                hold on;
                if strcmp(orientation, 'positive')
                    plot(locs_positive, modified_z_scored_signal(locs_positive), 'rv', 'MarkerFaceColor', 'r');
                    legend('LFP Signal', 'Detected Positive Peaks');
                else
                    plot(locs_negative, modified_z_scored_signal(locs_negative), 'bv', 'MarkerFaceColor', 'b');
                    legend('LFP Signal', 'Detected Negative Peaks');
                end
                title(['Detected Peaks - ', orientation, ' orientation (Signal #' num2str(k) ')']);
                xlabel('Samples');
                ylabel('Modified Z-scored LFP Signal');
                hold off;
                
                % Parameters for PQRST complex epoch extraction
                pre_peak_time = 200; % ms before peak, it can be modified
                post_peak_time = 400; % ms after peak, it can be modified
                pre_peak_samples = round(pre_peak_time * Fs / 1000);
                post_peak_samples = round(post_peak_time * Fs / 1000);
                epoch_length = pre_peak_samples + post_peak_samples + 1;
                
                % Extract PQRST epochs
                num_peaks = length(detected_peaks);
                pqrst_epochs = [];
                valid_epochs = 0;
                
                for i = 1:num_peaks
                    r_peak = detected_peaks(i);
                    start_idx = r_peak - pre_peak_samples;
                    end_idx = r_peak + post_peak_samples;
                    
                    % Check if epoch is within signal bounds
                    if start_idx >= 1 && end_idx <= length(LFP_signal)
                        epoch = LFP_signal(start_idx:end_idx);
                        valid_epochs = valid_epochs + 1;
                        pqrst_epochs(valid_epochs, :) = epoch;
                    end
                end
                
                % Only proceed if we have valid epochs
                if valid_epochs >= 3
                    % Average to create template
                    pqrst_template = mean(pqrst_epochs(1:valid_epochs, :), 1);
                    
                    % Plot PQRST template
                    figure;
                    plot(linspace(-pre_peak_time, post_peak_time, epoch_length), pqrst_template);
                    xlabel('Time (ms)');
                    ylabel('Amplitude');
                    title(['PQRST Complex Template (Signal #' num2str(k) ')']);
                    
                    % Subtract template from original signal
                    cleaned_LFP = LFP_signal;
                    
                    for i = 1:num_peaks
                        r_peak = detected_peaks(i);
                        start_idx = r_peak - pre_peak_samples;
                        end_idx = r_peak + post_peak_samples;
                        
                        if start_idx >= 1 && end_idx <= length(cleaned_LFP)
                            % Get the signal segment
                            signal_segment = cleaned_LFP(start_idx:end_idx);
                            
                            % Ensure template matches signal segment dimensions
                            if length(signal_segment) == length(pqrst_template)
                                % Make sure both are same orientation (row/column)
                                if size(signal_segment, 1) == 1 && size(pqrst_template, 1) > 1
                                    template_to_subtract = pqrst_template';
                                elseif size(signal_segment, 1) > 1 && size(pqrst_template, 1) == 1
                                    template_to_subtract = pqrst_template';
                                else
                                    template_to_subtract = pqrst_template;
                                end
                                
                                cleaned_LFP(start_idx:end_idx) = signal_segment - template_to_subtract;
                            else
                                fprintf('Warning: Template length mismatch for peak %d in signal %d\n', i, k);
                            end
                        end
                    end
                    
                    % Store cleaned LFP with triggers (overwrites the initialized copy)
                    d.cleaned_LFP{k} = [cleaned_LFP; d.butter_filtered{k}(2,:)];
                    
                    % Plot comparison
                    figure;
                    ax1 = subplot(2,1,1);
                    plot(LFP_signal);
                    title(['Original LFP Signal (Signal #' num2str(k) ')']);
                    xlabel('Samples');
                    ylabel('Amplitude');
                    
                    ax2 = subplot(2,1,2);
                    plot(cleaned_LFP);
                    title(['Cleaned LFP Signal (Signal #' num2str(k) ')']);
                    xlabel('Samples');
                    ylabel('Amplitude');
                    
                    linkaxes([ax1 ax2],'x');
                else
                    % Not enough valid epochs, keep butter_filtered version
                    fprintf('Warning: Not enough valid epochs for signal %d, keeping butter_filtered version\n', k);
                end
            else
                % Not enough peaks detected, keep butter_filtered version
                fprintf('Warning: Not enough peaks detected for signal %d, keeping butter_filtered version\n', k);
            end
        end
    end
    
    fprintf('Cardiac cleaning completed for specified condition.\n');
else
    fprintf('No cardiac cleaning applied - using butter_filtered data for all conditions.\n');
end

%% 4) Smoothing (optional step)

% Choose smoothing method: 'movmean' or 'sgolay'
smooth_method = 'movmean';  % options: 'movmean' or 'sgolay'

% Define the smoothing window size (in samples)
smooth_window = round(0.5 * Fs); % 0.5 needs to be customize depending on the filter chosen and your needs
fprintf('Applying smoothing (window = %.1f seconds) using %s...\n', smooth_window/Fs, smooth_method);

% If using Savitzky-Golay, define polynomial order
sgolay_order = 3;  % typical: 2-4, must be less than window size

for k = 1:4
    if ~isempty(d.cleaned_LFP{k})
        switch smooth_method
            case 'movmean'
                % Smooth the LFP signals using a moving average
                smoothed_lfp = movmean(d.cleaned_LFP{k}(1,:), smooth_window);

            case 'sgolay'
                % Ensure window length is odd and larger than polynomial order
                if mod(smooth_window,2) == 0
                    smooth_window = smooth_window + 1;
                end
                if smooth_window <= sgolay_order
                    error('SG window must be larger than polynomial order.');
                end
                % Apply Savitzky-Golay smoothing
                smoothed_lfp = sgolayfilt(d.cleaned_LFP{k}(1,:), sgolay_order, smooth_window);
        end

        % Keep original triggers
        d.smoothed_final{k} = [smoothed_lfp; d.raw{k}.triggers];
    else
        d.smoothed_final{k} = [];
    end
end


fprintf('LFP filtering completed.\n');

