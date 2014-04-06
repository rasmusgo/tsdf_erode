function plot_normpdf( mu, sigma, prob_range, xlims, ylims )
%% Plot probability density function as a histogram
cla; hold on

% Compute histogram
prob_range_limits = prob_range(1:end-1) + diff(prob_range)*0.5;
prob_cdf_values = [0, normcdf(prob_range_limits, mu, sigma), 1];
prob = diff(prob_cdf_values);

prob_scaled = prob./[1 diff(prob_range_limits) 1];
bar(prob_range, prob_scaled, 'FaceColor', [.5 .5 .5]);

%% Display expected value as a circle on the x-axis
plot(mu, 0, 'ob', 'LineWidth', 2, 'MarkerSize', 10);
Z = prob_range(:)' * prob(:) / sum(prob(:));
plot(Z, 0, 'xr', 'LineWidth', 2, 'MarkerSize', 10);

p = cumsum(prob(:));
i = find(p >= 0.5, 1, 'first');
if i > 1 && i < numel(p)
    y1 = p(i-1);
    y2 = p(i);
    x1 = prob_range_limits(i-1);
    x2 = prob_range_limits(i);
    plot(prob_range_limits, p(1:end-1), '-y')
    plot([x1 x2], [y1 y2], '-g', 'LineWidth', 3)
%     y = kx + m
%     k = (y2-y1) / (x2-x1);
%     m = y1 -k*x1;
%     x = (y - m)/k;
%     x = (y - y1 + k*x1)/k;
%     x = (y - y1)/k + x1;
%     x = (y - y1)*(x2-x1)/(y2-y1) + x1;
%     y = 0.5
%     x = (0.5 - y1)*(x2-x1)/(y2-y1) + x1;
    Z2 = (0.5 - y1)*(x2-x1)/(y2-y1) + x1;
else
    Z2 = prob_range(i);
end

plot(Z2, 0, '+g', 'LineWidth', 2, 'MarkerSize', 10);

%% Plot normal distribution on top
X = linspace(xlims(1), xlims(2), 300);
plot(X, normpdf(X, mu, sigma), 'LineWidth', 2);

xlim(xlims)
ylim(ylims)
trunc_dist = prob_range(end);
max_height = max(prob_scaled);
if max_height > ylims(2)
    title(sprintf('mu = %.1f sigma = %.1f trunc dist = %.1f (maximum bar height = %.2f)', ...
        mu, sigma, trunc_dist, max_height))
else
    title(sprintf('mu = %.1f sigma = %.1f trunc dist = %.1f', ...
        mu, sigma, trunc_dist))
end

hold off
end
