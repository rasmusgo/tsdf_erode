function plot_tsdf( Z, W, trunc_dist, P )
%PLOT_TSDF Show TSDF values and plot contour
    clf;
    subplot(1,3,1)

    % avoid plotting contours next to zero weight
    Z_nan = Z;
    Z_nan(W==0) = NaN;

    % plot distance
    subplot(1,2,1)
    hold on
    imagesc(Z, [-trunc_dist, trunc_dist])
    colormap(hot(1024))
    colorbar

    % Ground truth
    if size(P,1) == 2
        plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-k', 'LineWidth', 3)
        plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-w', 'LineWidth', 1)
    end

    % Indicate truncation distance
    contour(Z_nan, [-trunc_dist trunc_dist], '-g')

    % Plot zero crossing
    contour(Z_nan, [0 0], '-b', 'LineWidth', 2)

    axis equal xy tight
    hold off

    % plot weights
    subplot(1,2,2)
    hold on
    imagesc(W, [0,max(1, max(W(:)))])

    % Ground truth
    if size(P,1) == 2
        plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-k', 'LineWidth', 3)
        plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-w', 'LineWidth', 1)
    end

    % Indicate truncation distance
    contour(Z_nan, [-trunc_dist trunc_dist], '-g')

    % Plot zero crossing
    contour(Z_nan, [0 0], '-b', 'LineWidth', 2)
    colorbar

    axis equal xy tight
    hold off

    % set first subplot as active
    subplot(1,2,1)
end
