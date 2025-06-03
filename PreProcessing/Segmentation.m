%% LFP TRIGGER-BASED SEGMENTATION
% This script segments d.smoothed_final based on triggers
% Creates 5-second epochs after each trigger (to be modified depending on
% your experiment)

% Parameters
epoch_duration = 5; % seconds
Fs = data.BrainSenseTimeDomain(1).SampleRateInHz; % sampling rate
epoch_samples = round(epoch_duration * Fs); % samples per epoch

fprintf('Segmenting LFP data...\n');
fprintf('Epoch duration: %.1f seconds (%d samples)\n', epoch_duration, epoch_samples);
fprintf('Sampling rate: %d Hz\n', Fs);

% Initialize segmented data structure
segmented_data = struct();

% Process each channel
for k = 1:4
    if ~isempty(d.smoothed_final{k})
        fprintf('\nProcessing Channel %d...\n', k);
        
        % Extract LFP signal and triggers
        lfp_signal = d.smoothed_final{k}(1,:);
        trigger_signal = d.smoothed_final{k}(2,:);
        
        % Find trigger onsets (assuming triggers are binary: 0 = off, >0 = on)
        % Detect rising edges (transitions from 0 to >0)
        trigger_diff = diff([0, trigger_signal > 0]);
        trigger_onsets = find(trigger_diff == 1);
        
        % Alternative trigger detection methods (uncomment if needed):
        % Method 2: Find all non-zero trigger values
        % trigger_onsets = find(trigger_signal > 0);
        
        % Method 3: Find specific trigger values
        % trigger_onsets = find(trigger_signal == specific_trigger_value);
        
        fprintf('Found %d trigger onsets\n', length(trigger_onsets));
        
        if isempty(trigger_onsets)
            fprintf('Warning: No triggers found in channel %d\n', k);
            segmented_data.(['channel_' num2str(k)]) = [];
            continue;
        end
        
        % Initialize arrays for this channel
        valid_epochs = [];
        epoch_info = [];
        valid_count = 0;
        
        % Extract epochs
        for t = 1:length(trigger_onsets)
            trigger_time = trigger_onsets(t);
            epoch_start = trigger_time;
            epoch_end = trigger_time + epoch_samples - 1;
            
            % Check if epoch is within signal bounds
            if epoch_end <= length(lfp_signal)
                valid_count = valid_count + 1;
                
                % Extract epoch
                epoch_data = lfp_signal(epoch_start:epoch_end);
                valid_epochs(valid_count, :) = epoch_data;
                
                % Store epoch information
                epoch_info(valid_count).epoch_number = valid_count;
                epoch_info(valid_count).trigger_sample = trigger_time;
                epoch_info(valid_count).trigger_time_sec = (trigger_time - 1) / Fs;
                epoch_info(valid_count).trigger_value = trigger_signal(trigger_time);
                epoch_info(valid_count).epoch_start_sample = epoch_start;
                epoch_info(valid_count).epoch_end_sample = epoch_end;
                
                fprintf('Epoch %d: Trigger at sample %d (%.2f sec), value = %.1f\n', ...
                    valid_count, trigger_time, (trigger_time-1)/Fs, trigger_signal(trigger_time));
            else
                fprintf('Warning: Epoch %d extends beyond signal end (trigger at sample %d)\n', ...
                    t, trigger_time);
            end
        end
        
        % Store results for this channel
        if valid_count > 0
            segmented_data.(['channel_' num2str(k)]).epochs = valid_epochs;
            segmented_data.(['channel_' num2str(k)]).epoch_info = epoch_info;
            segmented_data.(['channel_' num2str(k)]).num_epochs = valid_count;
            segmented_data.(['channel_' num2str(k)]).epoch_duration_sec = epoch_duration;
            segmented_data.(['channel_' num2str(k)]).epoch_samples = epoch_samples;
            segmented_data.(['channel_' num2str(k)]).sampling_rate = Fs;
            
            fprintf('Successfully extracted %d valid epochs from channel %d\n', valid_count, k);
        else
            segmented_data.(['channel_' num2str(k)]) = [];
            fprintf('No valid epochs extracted from channel %d\n', k);
        end
    else
        fprintf('Channel %d is empty, skipping...\n', k);
        segmented_data.(['channel_' num2str(k)]) = [];
    end
end

%% OPTIONAL: Group epochs by trigger value (if you have different conditions)
fprintf('\n=== GROUPING BY TRIGGER VALUE ===\n');

% Initialize grouped data structure
grouped_data = struct();

for k = 1:4
    if ~isempty(segmented_data.(['channel_' num2str(k)]))
        channel_data = segmented_data.(['channel_' num2str(k)]);
        
        % Get unique trigger values
        trigger_values = [channel_data.epoch_info.trigger_value];
        unique_triggers = unique(trigger_values);
        
        fprintf('Channel %d has trigger values: %s\n', k, num2str(unique_triggers));
        
        % Group epochs by trigger value
        for trigger_val = unique_triggers
            % Find epochs with this trigger value
            matching_epochs = trigger_values == trigger_val;
            condition_name = ['condition_' num2str(trigger_val)];
            
            % Store grouped data
            if ~isfield(grouped_data, ['channel_' num2str(k)])
                grouped_data.(['channel_' num2str(k)]) = struct();
            end
            
            grouped_data.(['channel_' num2str(k)]).(condition_name).epochs = ...
                channel_data.epochs(matching_epochs, :);
            grouped_data.(['channel_' num2str(k)]).(condition_name).epoch_info = ...
                channel_data.epoch_info(matching_epochs);
            grouped_data.(['channel_' num2str(k)]).(condition_name).num_epochs = ...
                sum(matching_epochs);
            
            fprintf('  Condition %d: %d epochs\n', trigger_val, sum(matching_epochs));
        end
    end
end

%% VISUALIZATION (Optional)
fprintf('\n=== CREATING VISUALIZATION ===\n');

% Plot first few epochs from each channel
for k = 1:4
    if ~isempty(segmented_data.(['channel_' num2str(k)]))
        figure('Name', ['Channel ' num2str(k) ' - Segmented Epochs']);
        
        channel_data = segmented_data.(['channel_' num2str(k)]);
        epochs_to_plot = min(5, channel_data.num_epochs); % Plot first 5 epochs max
        
        time_vector = (0:epoch_samples-1) / Fs; % Time vector in seconds
        
        for ep = 1:epochs_to_plot
            subplot(epochs_to_plot, 1, ep);
            plot(time_vector, channel_data.epochs(ep, :));
            title(['Epoch ' num2str(ep) ' - Trigger Value: ' ...
                num2str(channel_data.epoch_info(ep).trigger_value)]);
            xlabel('Time (seconds)');
            ylabel('LFP Amplitude');
            grid on;
        end
        
        sgtitle(['Channel ' num2str(k) ' - First ' num2str(epochs_to_plot) ' Epochs']);
    end
end

%% SUMMARY
fprintf('\n=== SEGMENTATION SUMMARY ===\n');
fprintf('Epoch duration: %.1f seconds (%d samples)\n', epoch_duration, epoch_samples);
fprintf('Sampling rate: %d Hz\n', Fs);

total_epochs = 0;
for k = 1:4
    if ~isempty(segmented_data.(['channel_' num2str(k)]))
        num_epochs = segmented_data.(['channel_' num2str(k)]).num_epochs;
        fprintf('Channel %d: %d epochs\n', k, num_epochs);
        total_epochs = total_epochs + num_epochs;
    else
        fprintf('Channel %d: No data\n', k);
    end
end
fprintf('Total epochs extracted: %d\n', total_epochs);

fprintf('\nSegmentation completed!\n');
fprintf('Data structures created:\n');
fprintf('  - segmented_data: epochs organized by channel\n');
fprintf('  - grouped_data: epochs organized by channel and condition\n');
