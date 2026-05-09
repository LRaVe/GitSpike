function [C, spike_times] = f_spike_synchro(Spike_train1, Spike_train2, t_min, t_max)
    % Given two spike trains, this function calculates the coincidence of spikes between them.
    % Uses a greedy matching algorithm to pair spikes based on minimum distance.
    % Returns C: coincidence array (1 if matched and within tau, 0 otherwise)
    %        spike_times: corresponding spike times for each coincidence value
    
    % Ensure inputs are column vectors
    Spike_train1 = Spike_train1(:);
    Spike_train2 = Spike_train2(:);
    % Slice spikes within time window
    spike_train1_sliced = Spike_train1(Spike_train1 >= t_min & Spike_train1 <= t_max);
    spike_train2_sliced = Spike_train2(Spike_train2 >= t_min & Spike_train2 <= t_max);
    % Initialize tracking matrices
    n1 = length(spike_train1_sliced);
    n2 = length(spike_train2_sliced);
    if n1 == 0 && n2 == 0
        C = [];
        spike_times = [];
        return;
    end
    % Go through the spikes in train1 and find the closest spike in train2, using the adaptive tau for coincidence detection 
    % With the use of f_interval, and f_in_interval functions to determine if spikes are coincident
    C = zeros(n1, 1);  % Initialize coincidence vector
    spike_times = spike_train1_sliced;  % Corresponding spike times for C
    for i = 1:n1
        spike1 = spike_train1_sliced(i);
        tau1 = f_interval(spike_train1_sliced, spike1, t_min, t_max);
        for j = 1:n2
            spike2 = spike_train2_sliced(j);
            tau2 = f_interval(spike_train2_sliced, spike2, t_min, t_max);
            if f_in_interval(spike1, spike2, tau1, tau2)
                C(i) = 1;  % Mark as coincident
                break;  % Move to next spike in train1 after finding a match
            else
                C(i) = 0;  % Not coincident
            end
        end
    end
end

function min_interval = f_interval(spike_train, spike, t_min, t_max)
    % Calculate tau for adaptive coincidence detection
    % tau = min(forward_ISI, backward_ISI) / 2
    % Using eq. (15) and (16) for SPIKE-Synchronization
    
    spike_index = find(abs(spike_train - spike) < 1e-10, 1);
    
    if isempty(spike_index)
        min_interval = 0;
        return;
    end
    
    % Distance to previous spike
    if spike_index > 1
        prev_dist = spike_train(spike_index) - spike_train(spike_index - 1);
    else
        % For first spike: use ISI to next spike (or time window edge)
        if spike_index < length(spike_train)
            prev_dist = spike_train(spike_index + 1) - spike_train(spike_index);
        else
            prev_dist = (t_max - t_min);  % Single spike fallback
        end
    end
    
    % Distance to next spike
    if spike_index < length(spike_train)
        next_dist = spike_train(spike_index + 1) - spike_train(spike_index);
    else
        % For last spike: use ISI to previous spike (or time window edge)
        if spike_index > 1
            next_dist = spike_train(spike_index) - spike_train(spike_index - 1);
        else
            next_dist = (t_max - t_min);  % Single spike fallback
        end
    end
    
    min_interval = min(prev_dist, next_dist) / 2;
end

function check = f_in_interval(spike1, spike2, tau1, tau2)
    % Adaptive coincidence detection: eq. (15) and (16)
    % tau_ij = min(tau_i, tau_j)
    % Coincident if |spike1 - spike2| < tau_ij
    tau_ij = min(tau1, tau2);
    distance = abs(spike1 - spike2);
    
    if distance < tau_ij
        check = 1;
    else
        check = 0;
    end
end
