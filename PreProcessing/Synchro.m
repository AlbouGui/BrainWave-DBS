% Author: Albert Guillemette BSc MSc, 02.06.2025
%Synchronize triggers encoded in a excel table into raw LFP data
%To prepare the excel file, see the instruction on the Github

clear all; close all
% Load JSON output (contains data.BrainSenseTimeDomain)
[filename, pathname] = uigetfile('*.mat', 'Select the _output.mat file');
if filename == 0
    error('No file selected');
end
load(fullfile(pathname, filename)); % loads "data"

% Missing sample detection for DBS LFP data 

fprintf('Checking for missing samples in LFP data using packet-based approach...\n');

% Initialize tracking variables
total_missing_samples = 0;
recordings_with_gaps = 0;
corrected_recordings = 0;

% Assume typical LFP sampling rate (adjust based on your system)
expected_fs = 250; % Hz - modify this based on your actual sampling rate

for recording_idx = 1:length(data.BrainSenseTimeDomain)
    try
        % Get current recording
        current_data = data.BrainSenseTimeDomain(recording_idx);
        
        % Check if required fields exist
        if ~isfield(current_data, 'TimeDomainData') || isempty(current_data.TimeDomainData)
            continue;
        end
        
        % Extract data similar to the original function
        LFP_data = current_data.TimeDomainData;
        
        % Ensure data is a matrix (samples x channels)
        if size(LFP_data, 1) < size(LFP_data, 2)
            LFP_data = LFP_data';
        end
        
        % Get packet information
        if isfield(current_data, 'PacketGenTime') && isfield(current_data, 'PacketSampleCount')
            TicksInS = current_data.PacketGenTime;
            GlobalPacketSizes = current_data.PacketSampleCount;
        else
            % If packet info not available, try to estimate
            fprintf('Recording %d: Missing packet timing info, skipping detailed analysis\n', recording_idx);
            continue;
        end
        
        % Skip if insufficient packet data
        if length(TicksInS) < 2 || length(GlobalPacketSizes) < 1
            continue;
        end
        
        % Create theoretical time vector (similar to Thenaisie approach)
        time_start = TicksInS(1);
        time_end = TicksInS(end) + (GlobalPacketSizes(end)-1)/expected_fs;
        time = time_start:1/expected_fs:time_end;
        time = round(time, 3); % Round to millisecond precision
        
        % Create logical vector indicating which samples have been received
        isReceived = zeros(size(time, 2), 1);
        nPackets = numel(GlobalPacketSizes);
        
        % Process each packet (following Thenaisie logic)
        for packetId = 1:nPackets
            % Find closest time point for this packet
            timeTicksDistance = abs(time - TicksInS(packetId));
            [~, packetIdx] = min(timeTicksDistance);
            
            % Handle potential overlaps (from original code)
            if packetIdx <= length(isReceived) && isReceived(packetIdx) == 1
                packetIdx = packetIdx + 1;
            end
            
            % Mark samples as received for this packet
            packet_end_idx = min(length(isReceived), packetIdx + GlobalPacketSizes(packetId) - 1);
            if packetIdx <= length(isReceived)
                isReceived(packetIdx:packet_end_idx) = isReceived(packetIdx:packet_end_idx) + 1;
            end
        end
        
        % Count missing samples
        missing_samples = sum(isReceived == 0);
        if missing_samples > 0
            recordings_with_gaps = recordings_with_gaps + 1;
            total_missing_samples = total_missing_samples + missing_samples;
            
            fprintf('Recording %d: %d missing samples out of %d total (%.2f%%)\n', ...
                    recording_idx, missing_samples, length(isReceived), ...
                    100 * missing_samples / length(isReceived));
            
            % Optional: Create corrected data with NaNs (similar to original function)
            if missing_samples > 0
                % Create NaN-filled data matrix
                nChannels = size(LFP_data, 2);
                corrected_data = NaN(length(time), nChannels);
                
                % Fill in received samples
                received_mask = logical(isReceived);
                if sum(received_mask) == size(LFP_data, 1)
                    corrected_data(received_mask, :) = LFP_data;
                    
                    % Update the data structure
                    data.BrainSenseTimeDomain(recording_idx).TimeDomainData = corrected_data;
                    data.BrainSenseTimeDomain(recording_idx).CorrectedTime = time;
                    data.BrainSenseTimeDomain(recording_idx).MissingSamples = missing_samples;
                    
                    corrected_recordings = corrected_recordings + 1;
                    
                    fprintf('  → Corrected data created with NaN padding\n');
                else
                    fprintf('  → Warning: Sample count mismatch, cannot correct automatically\n');
                end
            end
            
            % Create visualization for first few problematic recordings
            if recordings_with_gaps <= 3
                figure('Position', [100 + 50*recordings_with_gaps, 100 + 50*recordings_with_gaps, 800, 400]);
                
                subplot(2,1,1);
                plot(isReceived, '.');
                yticks([0 1 2]);
                yticklabels({'not received', 'received', 'duplicate'});
                ylim([-0.5 2.5]);
                title(sprintf('Recording %d: Sample Reception Pattern', recording_idx));
                xlabel('Sample Index');
                ylabel('Reception Status');
                grid on;
                
                subplot(2,1,2);
                % Show gaps as red regions
                gap_indices = find(isReceived == 0);
                if ~isempty(gap_indices)
                    % Create binary plot showing gaps
                    gap_plot = ones(size(isReceived));
                    gap_plot(gap_indices) = 0;
                    area(gap_plot, 'FaceColor', [0.8 0.8 0.8]);
                    hold on;
                    area(1-gap_plot, 'FaceColor', [1 0.6 0.6]);
                    ylim([0 1]);
                    title(sprintf('Missing Data Regions (%d gaps)', length(gap_indices)));
                    xlabel('Sample Index');
                    ylabel('Data Presence');
                    legend({'Received', 'Missing'}, 'Location', 'best');
                end
                
                sgtitle(sprintf('Missing Sample Analysis - Recording %d', recording_idx), ...
                        'FontSize', 12, 'FontWeight', 'bold');
            end
        end
        
        % Check for duplicate samples (following original logic)
        duplicate_samples = sum(isReceived > 1);
        if duplicate_samples > 0
            fprintf('Recording %d: %d duplicate sample positions detected\n', ...
                    recording_idx, duplicate_samples);
        end
        
    catch ME
        fprintf('Error analyzing recording %d: %s\n', recording_idx, ME.message);
    end
end

% Summary report (similar to original function output)
fprintf('\n=== Missing Sample Detection Summary (Thenaisie Method) ===\n');
fprintf('Total recordings analyzed: %d\n', length(data.BrainSenseTimeDomain));
fprintf('Recordings with detected gaps: %d\n', recordings_with_gaps);
fprintf('Recordings corrected with NaN padding: %d\n', corrected_recordings);
fprintf('Total estimated missing samples: %d\n', total_missing_samples);

if total_missing_samples > 0
    fprintf('Estimated missing data duration: %.2f seconds\n', total_missing_samples / expected_fs);
    
    % Calculate total samples across all recordings
    total_samples = 0;
    for i = 1:length(data.BrainSenseTimeDomain)
        if isfield(data.BrainSenseTimeDomain(i), 'TimeDomainData')
            total_samples = total_samples + size(data.BrainSenseTimeDomain(i).TimeDomainData, 1);
        end
    end
    
    if total_samples > 0
        fprintf('Missing data percentage: %.3f%%\n', 100 * total_missing_samples / total_samples);
    end
    
    % Issue warning like in original function
    warning('Some samples were lost during these recordings and have been replaced with NaNs');
else
    fprintf('No missing samples detected - data appears complete!\n');
end

fprintf('Missing sample detection complete. Proceeding with trigger synchronization...\n\n');

% Load Excel file (DBS_times.xlsx)
[excelfile, excelpath] = uigetfile('*.xlsx', 'Select the DBS_times.xlsx file');
if excelfile == 0
    error('No Excel file selected');
end
[~, ~, dERt] = xlsread(fullfile(excelpath, excelfile));

% Column indices from Excel
col_rowjson = find(strcmp(dERt(1,:), 'Corresponding_row_in_json'));
col_sample = find(strcmp(dERt(1,:), 'DBS_samples'));
col_stn = find(strcmp(dERt(1,:), 'STN_side'));
col_stim = find(strcmp(dERt(1,:), 'DBS_stim_state'));

% Verify all columns were found
if isempty(col_rowjson) || isempty(col_sample) || isempty(col_stn) || isempty(col_stim)
    error('Required columns not found in Excel file');
end

% Prepare condition labels and mapping
cond_map = containers.Map({'LeftOff','RightOff','LeftOn','RightOn'}, 1:4);
labels = {'LstnOff', 'RstnOff', 'LstnOn', 'RstnOn'};
d.label = labels;

% Initialize data structure - one signal per condition with all triggers
d.raw = cell(4,1);
for j = 1:4
    d.raw{j} = struct('LFP', [], 'triggers', [], 'initialized', false);
end

% Counters for each condition
trial_counts = zeros(1, 4);

% Statistics tracking
total_processed = 0;
total_skipped = 0;

fprintf('Processing DBS data...\n');

% Loop through data rows in Excel (skip header row)
for i = 2:size(dERt, 1)
    try
        row_json = dERt{i, col_rowjson};
        sample = dERt{i, col_sample};
        side = dERt{i, col_stn};
        stim = dERt{i, col_stim};
        
        % Skip rows with missing values or 'N/A'
        if isempty(row_json) || isempty(sample) || isempty(side) || isempty(stim) || ...
           strcmp(row_json, 'N/A') || strcmp(sample, 'N/A') || strcmp(side, 'N/A') || strcmp(stim, 'N/A')
            total_skipped = total_skipped + 1;
            continue;
        end
        
        % Validate row_json index
        if row_json < 1 || row_json > length(data.BrainSenseTimeDomain)
            fprintf('Invalid JSON row index %d (row %d)\n', row_json, i);
            total_skipped = total_skipped + 1;
            continue;
        end
        
        % Get LFP signal
        signal = data.BrainSenseTimeDomain(row_json).TimeDomainData;
        
        % Ensure signal is a row vector
        if size(signal, 1) > size(signal, 2)
            signal = signal';
        end
        
        % Validate sample index
        if sample < 1 || sample > length(signal)
            fprintf('Sample %d out of bounds for signal length %d (row %d)\n', ...
                    sample, length(signal), i);
            total_skipped = total_skipped + 1;
            continue;
        end
        
        % Construct condition key
        key = [side stim]; % e.g., 'LeftOff', 'RightOn'
        
        % Check if key exists in map
        if ~isKey(cond_map, key)
            fprintf('Unknown condition: %s (row %d)\n', key, i);
            total_skipped = total_skipped + 1;
            continue;
        end
        
        % Get condition index
        idx = cond_map(key);
        
        % Initialize LFP and triggers for this condition if first trial
        if ~d.raw{idx}.initialized
            d.raw{idx}.LFP = signal;
            d.raw{idx}.triggers = zeros(1, length(signal));
            d.raw{idx}.initialized = true;
        end
        
        % Add trigger to the continuous LFP at the specified sample
        d.raw{idx}.triggers(sample) = 2;
        
        trial_counts(idx) = trial_counts(idx) + 1;
        total_processed = total_processed + 1;
        
        % Progress update every 10 trials
        if mod(total_processed, 10) == 0
            fprintf('Processed %d trials...\n', total_processed);
        end
        
    catch ME
        fprintf('Error processing row %d: %s\n', i, ME.message);
        total_skipped = total_skipped + 1;
    end
end

% Display summary
fprintf('\n=== Processing Complete ===\n');
fprintf('Total trials processed: %d\n', total_processed);
fprintf('Total trials skipped: %d\n', total_skipped);
fprintf('\nTrials per condition:\n');
for j = 1:4
    fprintf('  %s: %d triggers\n', labels{j}, trial_counts(j));
    if d.raw{j}.initialized
        trigger_positions = find(d.raw{j}.triggers > 0);
        fprintf('    Triggers at samples: [%s]\n', num2str(trigger_positions));
    end
end

% Now access your data simply as:
% d.raw{1}.LFP        - LFP signal for condition 1
% d.raw{1}.triggers   - All triggers for condition 1

% Plot LFP data with triggers for all conditions

% Create figure with subplots for each condition
figure('Position', [100, 100, 1200, 800]);

% Define colors for each condition
colors = {'b', 'r', 'g', 'm'}; % Blue, Red, Green, Magenta
condition_names = {'Left STN Off', 'Right STN Off', 'Left STN On', 'Right STN On'};

% Assuming sampling rate (adjust if you know the actual sampling rate)
fs = 250; % Hz - common for LFP data, adjust as needed

for cond = 1:4
    subplot(4, 1, cond);
    
    if d.raw{cond}.initialized && ~isempty(d.raw{cond}.LFP)
        % Get data
        lfp_signal = d.raw{cond}.LFP;
        triggers = d.raw{cond}.triggers;
        
        % Create time vector
        time = (0:length(lfp_signal)-1) / fs;
        
        % Plot LFP signal
        plot(time, lfp_signal, colors{cond}, 'LineWidth', 1);
        hold on;
        
        % Find and plot trigger positions
        trigger_indices = find(triggers > 0);
        if ~isempty(trigger_indices)
            trigger_times = trigger_indices / fs;
            trigger_values = lfp_signal(trigger_indices);
            
            % Plot triggers as red markers
            plot(trigger_times, trigger_values, 'ro', 'MarkerSize', 8, ...
                 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black', 'LineWidth', 2);
            
            % Add vertical lines at trigger positions
            for i = 1:length(trigger_times)
                line([trigger_times(i), trigger_times(i)], ylim, ...
                     'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
            end
        end
        
        % Formatting
        title(sprintf('%s (%s) - %d triggers', condition_names{cond}, labels{cond}, length(trigger_indices)));
        xlabel('Time (seconds)');
        ylabel('LFP Amplitude (µV)');
        grid on;
        
        % Add legend only for first subplot
        if cond == 1
            legend({'LFP Signal', 'Triggers'}, 'Location', 'best');
        end
        
        hold off;
    else
        % Empty condition
        text(0.5, 0.5, sprintf('No data for %s', condition_names{cond}), ...
             'HorizontalAlignment', 'center', 'FontSize', 12);
        title(sprintf('%s (%s) - No data', condition_names{cond}, labels{cond}));
    end
end

% Add overall title
sgtitle('LFP Signals with DBS Triggers by Condition', 'FontSize', 16, 'FontWeight', 'bold');

% Adjust spacing between subplots
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.1, 0.1, 0.8, 0.8]);

fprintf('\nPlot created! Red circles and dashed lines indicate trigger positions.\n');

% Optional: Create a zoomed-in view around triggers
figure('Position', [150, 150, 1000, 600]);

for cond = 1:4
    if d.raw{cond}.initialized && ~isempty(d.raw{cond}.LFP)
        triggers = d.raw{cond}.triggers;
        trigger_indices = find(triggers > 0);
        
        if ~isempty(trigger_indices)
            subplot(2, 2, cond);
            
            % Take first trigger for zoomed view (adjust window as needed)
            center_sample = trigger_indices(1);
            window_samples = round(2 * fs); % 2 seconds around trigger
            
            start_idx = max(1, center_sample - window_samples);
            end_idx = min(length(d.raw{cond}.LFP), center_sample + window_samples);
            
            % Extract segment
            time_segment = (start_idx:end_idx) / fs;
            lfp_segment = d.raw{cond}.LFP(start_idx:end_idx);
            
            % Plot segment
            plot(time_segment, lfp_segment, colors{cond}, 'LineWidth', 1.5);
            axis off
            hold on;
            
            % Mark trigger
            trigger_time = center_sample / fs;
            trigger_value = d.raw{cond}.LFP(center_sample);
            %plot(trigger_time, trigger_value, 'ro', 'MarkerSize', 10, ...
                % 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black', 'LineWidth', 2);
            
            % Add vertical line at trigger
            line([trigger_time, trigger_time], ylim, ...
                 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
            
            title(sprintf('%s - Zoomed View', labels{cond}));
            xlabel('Time (seconds)');
            ylabel('LFP Amplitude (µV)');
            grid off;
            hold off;
        end
    end
end

sgtitle('Zoomed View Around First Trigger (±2 seconds)', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('Zoomed plots created showing detail around trigger events.\n');


