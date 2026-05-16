function spike_train_order()
    % Compute the spike-train-order value for a set of spike trains
    % This function loads spike data, computes pairwise train orders, 
    % and outputs the overall spike-train-order value.

    % Ensure the whole repository is on the MATLAB path
    thisFile = mfilename('fullpath');
    if ~isempty(thisFile)
        repoRoot = fileparts(fileparts(thisFile));
        addpath(genpath(repoRoot));
    else
        addpath(genpath('..'));
    end

    % Define time window and initialize spike train data
    tmin=0;
    tmax=10;
    num_trains=4;
    spikes=cell(1,num_trains);
    %spikes{1} = [0.0001 0.7142];
    %spikes{2} = [0.2858 0.9999];                   
    %spikes{3} = [0.1429 0.8571];
    spikes{1}=[0 1.9 3.9 7 10];
    spikes{2}=[0 2 7.1 9 10];
    spikes{3}=[0 2.1 4.1 6.9 10];
    spikes{4}=[0 2.2 6.8 7.1 10];
    number_spikes=sum(cellfun(@length, spikes));
    
    % Add auxiliary spikes at boundaries (tmin and tmax)
    [spikes, ~, ~] = add_auxiliary_spikes(spikes, tmin, tmax);
    
    % Compute pairwise spike-train-order relationships
    [results,order_matrix]=order_trains(tmin,tmax,spikes);
    
    % Calculate the overall spike-train-order value F
    F = compute_spike_train_order_value(spikes, results, number_spikes);

    if nargout==0
        spike_train_order_profile=plot_spike_train_order(spikes,results,order_matrix,F,tmin,tmax);
    end

    % Display results to console
    fprintf('Spike train order F = %.4f\n', F);
    disp('Spike train order profile: ');
    disp(spike_train_order_profile);
    disp('Spike train order matrix:');
    disp(order_matrix);
end