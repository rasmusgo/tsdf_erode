function plot_pdf( sample, prob_range )
%% Plot probability density function as a histogram
cla; hold on

h = bar(prob_range, sample(:));

% Change color of the bars to match the corresponding distance
ch = get(h,'Children');
fvd = get(ch,'Faces');
fvcd = get(ch,'FaceVertexCData');
for i = 1:numel(prob_range)
    fvcd(fvd(i,:)) = i;
end
set(ch,'FaceVertexCData',fvcd)

colormap hot
ylim([-0.2, 1.2])

%% Display expected value as a vertical blue line
Z = prob_range(:)' * sample(:) / sum(sample(:));
plot([Z,Z], [-0.1, 1.1], '-b', 'LineWidth', 3);

hold off
end
