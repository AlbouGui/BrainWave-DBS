%% VISUALIZATION
par = {'raw', 'no_60Hz', 'butter_filtered', 'smoothed_final'};

fprintf('Creating visualization plots...\n');

for pn = 1:numel(par)
    figure('units','normalized','outerposition',[0 0 1 1],'color','w','visible','on');
    
    for p = 1:4 % 4 conditions
        if strcmp(par{pn}, 'raw')
            % Handle raw data structure
            if d.raw{p}.initialized && ~isempty(d.raw{p}.LFP)
                signal_data = d.raw{p}.LFP;
                trigger_data = d.raw{p}.triggers;
            else
                signal_data = [];
                trigger_data = [];
            end
        else
            % Handle filtered data structures
            if ~isempty(d.(par{pn}){p})
                signal_data = d.(par{pn}){p}(1,:);
                trigger_data = d.(par{pn}){p}(2,:);
            else
                signal_data = [];
                trigger_data = [];
            end
        end
        
        subplot(4,1,p)
        hold on
        
        if ~isempty(signal_data)
            % Create time vector in minutes
            time_min = (1:length(signal_data))/Fs/60;
            
            % Plot LFP signal
            plot(time_min, signal_data, 'color', 'b', 'LineWidth', 1)
            
            % Find and plot trigger positions
            temp_mark = find(trigger_data == 2);
            if ~isempty(temp_mark)
                plot(temp_mark/Fs/60, repmat(max(signal_data)*0.9, length(temp_mark), 1), 'r*', 'MarkerSize', 8)
            end
            
            % Formatting
            xlabel('Time (min)')
            ylabel('LFP Amplitude')
            title(labels{p}, 'interpreter', 'none')
            grid on
            
            % Add some statistics to title
            if ~isempty(temp_mark)
                title(sprintf('%s (%d triggers)', labels{p}, length(temp_mark)), 'interpreter', 'none')
            end
        else
            % Empty condition
            text(0.5, 0.5, sprintf('No data for %s', labels{p}), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'Units', 'normalized');
            title(labels{p}, 'interpreter', 'none')
        end
        
        hold off
    end
    
    % Add overall title
    sgtitle(par{pn}, 'interpreter', 'none', 'FontSize', 16, 'FontWeight', 'bold')
    
    % Save figure option
    saveas(gcf, sprintf('LFP_%s_filtering.png', par{pn}));
end

fprintf('Filtering complete! Created %d visualization figures.\n', length(par));

% Display filtering summary
fprintf('\n=== Filtering Summary ===\n');
for k = 1:4
    if d.raw{k}.initialized
        fprintf('%s:\n', labels{k});
        fprintf('  Signal length: %d samples (%.1f minutes)\n', ...
                length(d.raw{k}.LFP), length(d.raw{k}.LFP)/Fs/60);
        fprintf('  Number of triggers: %d\n', sum(d.raw{k}.triggers > 0));
    else
        fprintf('%s: No data\n', labels{k});
    end
end
