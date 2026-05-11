function plot_spike_train_order(spikes,orders,order_matrix,F,tmin,tmax)
    % Visualize spike-train-order results with two plots:
    % 1) Temporal ordering of all spikes with the F metric
    % 2) Heatmap of pairwise train ordering relationships
    
    n=length(spikes);
    time = horzcat(spikes{:});
    value = horzcat(orders{:});

    [sortedTimes,orderInd]=sort(time);
    sortedOrders=value(orderInd);

    % Create first figure: Temporal spike ordering
    figure;
    hold on;
    grid on;

    plot([tmin,tmax],[F,F], '-', 'Color', 'red', 'LineWidth', 1);
    plot(sortedTimes,sortedOrders,'-o','Color', 'blue', 'LineWidth', 1.5, 'MarkerSize', 6);

    xlim([tmin,tmax]);
    ylim([-1.1,1.1]);
    title(sprintf('Spike train order F = %.4g', F));
    yticks([-1,0,1]);
    hold off;

    % Create second figure: Pairwise ordering heatmap
    figure;
    hold on;
    matrix_min = min(order_matrix(:));
    matrix_max = max(order_matrix(:));
    imagesc(order_matrix, [matrix_min matrix_max]);
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
    title(sprintf('Pairwise train order matrix F = %g', F));
    hold off;
end