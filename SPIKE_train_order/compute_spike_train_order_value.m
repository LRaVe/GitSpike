function F = compute_spike_train_order_value(spikes, orders, number_spikes)
    % Calculate the overall spike-train-order value F
    
    % Concatenate all spike times and their corresponding order values
    time = horzcat(spikes{:});
    value = horzcat(orders{:});

    [~,orderInd]=sort(time);
    sortedOrders=value(orderInd);
    F=sum(sortedOrders);

    % Normalize by number of spikes
    if number_spikes~=0
        F=F/number_spikes;
    end
end