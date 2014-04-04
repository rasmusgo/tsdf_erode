function plot_open_closed( Zopen, Zclosed, trunc_dist, P )

%% Render for debugging
clf;

% plot open hypothesis
subplot(1,2,1)
title 'Open hypothesis'
hold on
imagesc(Zopen, [-trunc_dist, trunc_dist])
colormap(hot(1024))
colorbar

% Ground truth polygon
if size(P,1) == 2
    plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-k', 'LineWidth', 3)
    plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-w', 'LineWidth', 1)
end

% Indicate truncation distance
contour(Zopen, [-trunc_dist trunc_dist], '-g')

axis equal xy tight
hold off

% plot closed hypothesis
subplot(1,2,2)
title 'Closed hypothesis'
hold on
imagesc(Zclosed, [-trunc_dist trunc_dist])

% Ground truth polygon
if size(P,1) == 2
    plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-k', 'LineWidth', 3)
    plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-w', 'LineWidth', 1)
end

% Indicate truncation distance
contour(Zclosed, [-trunc_dist trunc_dist], '-g')

% Plot zero crossing
contour(Zclosed, [0 0], '-b', 'LineWidth', 2)
colorbar

axis equal xy tight
hold off
end
