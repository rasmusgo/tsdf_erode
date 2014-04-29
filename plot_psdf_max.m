function plot_psdf_max( prob, prob_range, P )
%% Animate probability signed distance field
cla

Np = size(prob,3);
% for i = 1:Np
%     imagesc(prob(:,:,i)*i, [0 Np])
%     axis equal xy tight
%     colorbar
%     drawnow
% end

%% Compute values with highest probability
[m, ind] = max(prob, [], 3);
Z = reshape(prob_range(reshape(ind, numel(ind), 1)), size(ind));

imagesc(Z, [prob_range(1) prob_range(end)])
hold on
colorbar
colormap hot

% Ground truth
if size(P,1) == 2
    plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-k', 'LineWidth', 3)
    plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-w', 'LineWidth', 1)
end

% Indicate truncation distance
contour(Z, [prob_range(1) prob_range(end)], '-g')

% Plot zero crossing
contour(Z, [0 0], '-b', 'LineWidth', 2)

axis equal xy tight
hold off

end
