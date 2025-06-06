% Author: Albert Guillemette BSc MSc, 02.06.2025
% VISUALIZATION (LFP & SPECTROGRAMS) with triggers
par = {'raw', 'no_60Hz', 'butter_filtered', 'smoothed_final'};
fprintf('Creating visualization plots with spectrograms...\n');

% Spectrogram parameters (can be modified)
window_length = 2 * Fs;  % 2 second window
overlap = 0.75;          % 75% overlap
nfft = 2^nextpow2(window_length);  % FFT points
max_freq = 100;          % Maximum frequency to display (Hz)

for pn = 1:numel(par)
    % Create figure with more space for spectrograms
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
        
        if ~isempty(signal_data)
            % TIME DOMAIN PLOT (left column)
            subplot(4,2,2*p-1)
            hold on
            
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
            if ~isempty(temp_mark)
                title(sprintf('%s - Time Domain (%d triggers)', labels{p}, length(temp_mark)), 'interpreter', 'none')
            else
                title(sprintf('%s - Time Domain', labels{p}), 'interpreter', 'none')
            end
            grid on
            hold off
            
            % SPECTROGRAM PLOT (right column)
            subplot(4,2,2*p)
            
            % Calculate spectrogram
            [S, F, T, P] = spectrogram(signal_data, window_length, round(window_length*overlap), nfft, Fs);
            
            % Convert to time in minutes
            T_min = T/60;
            
            % Limit frequency range
            freq_idx = F <= max_freq;
            F_display = F(freq_idx);
            P_display = P(freq_idx, :);
            
            % Plot spectrogram
            imagesc(T_min, F_display, 10*log10(P_display));
            axis xy;
            colormap(parula);
            colorbar;
            
            % Add trigger lines if they exist
            if ~isempty(temp_mark)
                hold on
                trigger_times_min = temp_mark/Fs/60;
                for t = 1:length(trigger_times_min)
                    line([trigger_times_min(t) trigger_times_min(t)], [0 max_freq], ...
                         'Color', 'white', 'LineWidth', 2, 'LineStyle', '--');
                end
                hold off
            end
            
            xlabel('Time (min)')
            ylabel('Frequency (Hz)')
            title(sprintf('%s - Spectrogram', labels{p}), 'interpreter', 'none')
            caxis([-40 10]); % Adjust power scale as needed
            
        else
            % Empty condition - both subplots
            subplot(4,2,2*p-1)
            text(0.5, 0.5, sprintf('No data for %s', labels{p}), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'Units', 'normalized');
            title(sprintf('%s - Time Domain', labels{p}), 'interpreter', 'none')
            
            subplot(4,2,2*p)
            text(0.5, 0.5, sprintf('No data for %s', labels{p}), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'Units', 'normalized');
            title(sprintf('%s - Spectrogram', labels{p}), 'interpreter', 'none')
        end
    end
    
    % Add overall title
    sgtitle(sprintf('%s - Time Domain & Spectrograms', par{pn}), 'interpreter', 'none', 'FontSize', 16, 'FontWeight', 'bold')
    
    % Save figure
    saveas(gcf, sprintf('LFP_%s_with_spectrograms.png', par{pn}));
end

%% ADDITIONAL: Create separate detailed spectrogram figures
fprintf('Creating detailed spectrogram figures...\n');

for pn = 1:numel(par)
    figure('units','normalized','outerposition',[0 0 1 1],'color','w','visible','on');
    
    for p = 1:4
        if strcmp(par{pn}, 'raw')
            if d.raw{p}.initialized && ~isempty(d.raw{p}.LFP)
                signal_data = d.raw{p}.LFP;
                trigger_data = d.raw{p}.triggers;
            else
                signal_data = [];
                trigger_data = [];
            end
        else
            if ~isempty(d.(par{pn}){p})
                signal_data = d.(par{pn}){p}(1,:);
                trigger_data = d.(par{pn}){p}(2,:);
            else
                signal_data = [];
                trigger_data = [];
            end
        end
        
        subplot(2,2,p)
        
        if ~isempty(signal_data)
            % Calculate high-resolution spectrogram
            [S, F, T, P] = spectrogram(signal_data, window_length, round(window_length*overlap), nfft, Fs);
            
            % Convert to time in minutes
            T_min = T/60;
            
            % Limit frequency range for better visualization
            freq_idx = F <= max_freq;
            F_display = F(freq_idx);
            P_display = P(freq_idx, :);
            
            % Plot spectrogram with better resolution
            imagesc(T_min, F_display, 10*log10(P_display));
            axis xy;
            colormap(parula);
            c = colorbar;
            c.Label.String = 'Power (dB)';
            
            % Add trigger markers
            temp_mark = find(trigger_data == 2);
            if ~isempty(temp_mark)
                hold on
                trigger_times_min = temp_mark/Fs/60;
                for t = 1:length(trigger_times_min)
                    line([trigger_times_min(t) trigger_times_min(t)], [0 max_freq], ...
                         'Color', 'white', 'LineWidth', 2, 'LineStyle', '--');
                end
                % Add trigger count to title
                title(sprintf('%s (%d triggers)', labels{p}, length(temp_mark)), 'interpreter', 'none')
                hold off
            else
                title(labels{p}, 'interpreter', 'none')
            end
            
            xlabel('Time (min)')
            ylabel('Frequency (Hz)')
            
            % Set consistent color scale
            caxis([-40 0]);
            
        else
            text(0.5, 0.5, sprintf('No data for %s', labels{p}), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, 'Units', 'normalized');
            title(labels{p}, 'interpreter', 'none')
        end
    end
    
    sgtitle(sprintf('%s - Detailed Spectrograms', par{pn}), 'interpreter', 'none', 'FontSize', 16, 'FontWeight', 'bold')
    saveas(gcf, sprintf('LFP_%s_spectrograms_detailed.png', par{pn}));
end

fprintf('Filtering and spectrogram visualization complete! Created %d combined and %d detailed spectrogram figures.\n', length(par), length(par));

% Display filtering summary
fprintf('\n=== Filtering Summary ===\n');
for k = 1:4
    if d.raw{k}.initialized
        fprintf('%s:\n', labels{k});
        fprintf('  Signal length: %d samples (%.1f minutes)\n', ...
                length(d.raw{k}.LFP), length(d.raw{k}.LFP)/Fs/60);
        fprintf('  Number of triggers: %d\n', sum(d.raw{k}.triggers > 0));
        
        % Add spectral analysis summary
        if ~isempty(d.raw{k}.LFP)
            % Calculate power spectral density
            [pxx, f] = pwelch(d.raw{k}.LFP, window_length, round(window_length*overlap), nfft, Fs);
            
            % Find dominant frequencies
            [~, peak_idx] = findpeaks(pxx, 'MinPeakHeight', max(pxx)*0.1, 'NPeaks', 3);
            dominant_freqs = f(peak_idx);
            
            fprintf('  Dominant frequencies: ');
            if ~isempty(dominant_freqs)
                fprintf('%.1f', dominant_freqs(1));
                for i = 2:length(dominant_freqs)
                    fprintf(', %.1f', dominant_freqs(i));
                end
                fprintf(' Hz\n');
            else
                fprintf('None detected\n');
            end
        end
    else
        fprintf('%s: No data\n', labels{k});
    end
end

fprintf('All visualizations complete!\n');