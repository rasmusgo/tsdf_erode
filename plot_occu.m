function plot_occu( logodds, P )
%PLOT_OCCU Show logodds values and plot contour
clf;

% plot logodds
hold on
imagesc(logodds, [log(0.01), log(100)])
colormap(jet(1024))
colorbar

% Ground truth
if size(P,1) == 2
    plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-k', 'LineWidth', 3)
    plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-w', 'LineWidth', 1)
end

% Plot 50% probability crossing
contour(logodds, [log(1) log(1)], '-m', 'LineWidth', 2)

axis equal xy tight
hold off

end
