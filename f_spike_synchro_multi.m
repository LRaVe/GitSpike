function [C_matrix,C_global] = f_spike_synchro_multi(st, t_min, t_max)
    % Multivariate SPIKE-Synchronization following equations (17-19) from the paper
    n_trains = length(st);
    C_matrix = zeros(n_trains, n_trains);
    %Compute pairwise coincidence matrices for all pairs of spike trains, without matrix for now
    for i = 1:n_trains
        for j = [1:n_trains]  % Compare train i to all other trains
            if i == j
                C_matrix(i, j) = 1;  % Perfect synchronization with itself
                continue;  % Skip self-comparison
            end
            [C_ij, times_ij] = f_spike_synchro(st{i}, st{j}, t_min, t_max);
            %Update C_matrix
            for k = 1:length(C_ij)
                C_matrix(i, j) = C_matrix(i, j) + C_ij(k);
            end
            C_matrix(i, j) = C_matrix(i, j) / (length(times_ij));  % Average over spikes
        end
    end
    %calculate global SPIKE-Synchronization index C_global as the mean of the upper triangle of C_matrix (excluding diagonal)
    C_global = mean(C_matrix(triu(true(size(C_matrix)), 1)));
end