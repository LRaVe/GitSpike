function fig_handles = f_spike_synchro_plot(fig_input, data, params)
    % Plotting function for SPIKE-Synchronization results
    % This function can be used to visualize the results of the SPIKE-Synchronization analysis
    % It can be called after running f_spike_synchro_multi to display the pairwise
    % coincidence matrix and the global synchronization index
    % Extract data
    C_matrix = data.C_matrix;
    C_global = data.C_global;
    spike_synchro_data = data.spike_synchro_data; % Cell array with [time, C_value] pairs
    
    % Plot pairwise coincidence matrix
    figure;
    imagesc(C_matrix);
    colorbar;
    colormap(jet);
    xlabel('Spike Train Index');
    ylabel('Spike Train Index');
    title(sprintf('Pairwise Coincidence Matrix (C_{global} = %.4f)', C_global));

    % Profile plot of synchronization over time for each train pair
    figure;
    n_trains = size(spike_synchro_data, 1);
    for i = 1:n_trains
        for j = 1:n_trains
            if i ~= j
                data_pair = spike_synchro_data{i, j};
                if ~isempty(data_pair)
                    times = data_pair(:, 1);  % Time values
                    C_vals = data_pair(:, 2); % C values
                    plot(times, C_vals, 'o', 'DisplayName', sprintf('Train %d vs %d', i, j));
                    hold on;
                end
            end
        end
    end
    grid on;
    xlabel('Time (ms)');
    ylabel('Coincidence Value (C)');
    title('SPIKE-Synchronization Over Time');
    legend('show');
    hold off;

end 