function plot_spike_order(spikes,SO_matrix,tmin,tmax,threshold) 
    % Visualize spike-order results for all spikes across spike trains
    % Displays individual spike orderings and the overall D metric
    
    % Compute spike-order values
    n=length(spikes);
    orders=order_spikes(tmin,tmax,spikes);
    
    time = horzcat(spikes{:}); 
    value = horzcat(orders{:}); 

    [sortedTimes,orderInd]=sort(time);
    sortedOrders=value(orderInd);

    D=sum(sortedOrders);
    
    % Round near-zero values to zero for cleaner display
    if abs(D) < threshold
        D = 0;
    end

    % Create first figure: temporal spike ordering
    figure;
    hold on;
    grid on;

    plot([tmin,tmax],[D,D], '-', 'Color', 'red', 'LineWidth', 1);
    plot(sortedTimes,sortedOrders,'-o','Color', 'blue', 'LineWidth', 1.5, 'MarkerSize', 6);

    xlim([tmin,tmax]);
    ylim([-1.1,1.1]);
    title(sprintf('Spike order D = %.4g', D));
    yticks([-1,0,1]);
    hold off;

    % Create second figure: pairwise ordering matrix
    figure;
    hold on;
    
    imagesc(SO_matrix, [-1 1]);
    colormap(jet);
    colorbar;
    axis equal;
    xlim([0.5 n+0.5]);
    ylim([0.5 n+0.5]);
    set(gca, 'XDir', 'normal');
    set(gca, 'YDir', 'reverse');
    set(gca, 'XTick', 1:n, 'YTick', 1:n);
    xlabel('Train index');
    ylabel('Train index');
    title(sprintf('Spike order matrix D = %g', D));
    hold off;
end
