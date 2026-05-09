function plot = f_spike_synchro_plot(plot, data, params)
    % Plotting function for SPIKE-Synchronization results
    % This function can be used to visualize the results of the SPIKE-Synchronization analysis
    % It can be called after running f_spike_synchro_multi to display the pair
    %wise coincidence matrix and the global synchronization index
    % Extract data
    C_matrix = data.C_matrix;
    C_global = data.C_global;
    % Plot pairwise coincidence matrix
    imagesc(C_matrix);
    colorbar;
    colormap(jet);
    title('Pairwise Coincidence Matrix');
    xlabel('Spike Train Index');
    ylabel('Spike Train Index');
    title(sprintf('Pairwise Coincidence Matrix (C_{global} = %.4f)', C_global));
end 