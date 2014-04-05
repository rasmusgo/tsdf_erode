function plot_normpdf( mu, sigma, prob_range, xlims, ylims )
%% Plot probability density function as a histogram
cla; hold on

% Compute histogram
prob_range_limits = prob_range(1:end-1) + diff(prob_range)*0.5;
prob_cdf_values = [0, normcdf(prob_range_limits, mu, sigma), 1];
prob = diff(prob_cdf_values);

prob_scaled = prob./[1 diff(prob_range_limits) 1];
bar(prob_range, prob_scaled, 'FaceColor', [.5 .5 .5]);

%% Display expected value as a vertical blue line
plot(mu, 0, 'ob', 'LineWidth', 3);
Z = prob_range(:)' * prob(:) / sum(prob(:));
plot(Z, 0, 'xr', 'LineWidth', 3);

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
