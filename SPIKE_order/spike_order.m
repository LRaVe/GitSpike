function results=spike_order()
    % Compute and display spike-order analysis for test spike trains
    % This script demonstrates the spike-order metric calculation and visualization
    
    threshold=1e-10;

    % Ensure spike_common is on the MATLAB path
    thisFile = mfilename('fullpath');
    if ~isempty(thisFile)
        repoRoot = fileparts(fileparts(thisFile));
        addpath(fullfile(repoRoot, 'spike_common'));
    else
        addpath('../spike_common');
    end

    % Define observation window and initialize spike trains
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
    
    % Add auxiliary spikes at boundaries
    [spikes, ~, ~] = add_auxiliary_spikes(spikes, tmin, tmax);
    
    % Compute spike-order values
    [results, SO_matrix] = order_spikes(tmin,tmax,spikes);

    if nargout==0
        % Visualize results if no output is requested
        plot_spike_order(spikes,SO_matrix,tmin,tmax,threshold);
    end

    disp('Spike-order results:');
    disp(results);
end