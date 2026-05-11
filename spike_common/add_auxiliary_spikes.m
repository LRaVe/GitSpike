function [spikes, aux_begin, aux_end] = add_auxiliary_spikes(spikes, t_min, t_max)
    % Add auxiliary spike times at the beginning and end of the observation window

    % ===============================================
    % ====== Handle cell array of spike trains ======
    % ===============================================

    if iscell(spikes)
        aux_begin = zeros(1, length(spikes));
        aux_end = zeros(1, length(spikes));
        for i = 1:length(spikes)
            [spikes{i}, aux_begin(i), aux_end(i)] = process_single_train(spikes{i}, t_min, t_max);
        end
        return;
    end

    % =======================================
    % ====== Handle single spike train ======
    % =======================================

    [spikes, aux_begin, aux_end] = process_single_train(spikes, t_min, t_max);
end

function [train, aux_begin, aux_end] = process_single_train(train, t_min, t_max)
    aux_begin = 0;
    aux_end = 0;

    train = sort(unique(train));

    % Add auxiliary spike at the beginning if needed
    if train(1) > t_min

        if length(train) >= 2
            aux = train(1) - max(train(1)-t_min, train(2)-train(1));
        else
            aux = t_min;
        end

        train = [aux train];
        aux_begin = 1;
    end

    % Add auxiliary spike at the end if needed
    if train(end) < t_max

        if length(train) >= 2
            aux = train(end) + max(t_max-train(end), train(end)-train(end-1));
        else
            aux = t_max;
        end

        train = [train aux];
        aux_end = 1;
    end
end