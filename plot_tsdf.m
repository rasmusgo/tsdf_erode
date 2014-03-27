function plot_tsdf( Z, W, trunc_dist )
%PLOT_TSDF Show TSDF values and plot contour
    clf;
    hold on
    imagesc(Z, [-trunc_dist, trunc_dist])
    colormap(hot(1024))
    colorbar
    Z(W==0) = NaN;
    contour(Z, [-trunc_dist trunc_dist], '-g')
    contour(Z, [0 0], '-b', 'LineWidth', 2)
    %colormap hot
    axis equal xy tight
    hold off
end

