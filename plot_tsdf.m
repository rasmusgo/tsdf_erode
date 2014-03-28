function plot_tsdf( Z, W, trunc_dist )
%PLOT_TSDF Show TSDF values and plot contour
    clf;

    % avoid plotting next to zero weight
    Z(W==0) = NaN;

    % plot distance
    subplot(1,2,1)
    hold on
    imagesc(Z, [-trunc_dist, trunc_dist])
    colormap(hot(1024))
    colorbar
    contour(Z, [-trunc_dist trunc_dist], '-g')
    contour(Z, [0 0], '-b', 'LineWidth', 2)
    %colormap hot
    axis equal xy tight
    hold off

    % plot weights
    subplot(1,2,2)
    hold on
    imagesc(W, [0,max(1, max(W(:)))])
    contour(Z, [-trunc_dist trunc_dist], '-g')
    contour(Z, [0 0], '-b', 'LineWidth', 2)
    %colormap(gray)
    colorbar
    axis equal xy tight
    hold off

    % set first subplot as active
    subplot(1,2,1)
end

