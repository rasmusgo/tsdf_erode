function plot_psdf( prob, trunc_dist, P )
%% Animate probability signed distance field
clf

Np = size(prob,3);
% for i = 1:Np
%     imagesc(prob(:,:,i)*i, [0 Np])
%     axis equal xy tight
%     colorbar
%     drawnow
% end

%% Compute expected value
prob_range = linspace(-trunc_dist, trunc_dist, Np);
Z = sum(bsxfun(@times, reshape(prob_range, 1,1,Np), prob), 3) ...
    ./ sum(prob, 3);

imagesc(Z, [-trunc_dist trunc_dist])
hold on
colorbar
colormap hot

% Ground truth
if size(P,1) == 2
    plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-k', 'LineWidth', 3)
    plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-w', 'LineWidth', 1)
end

% Indicate truncation distance
contour(Z, [-trunc_dist trunc_dist], '-g')

% Plot zero crossing
contour(Z, [0 0], '-b', 'LineWidth', 2)

axis equal xy tight
hold off

end
