%% Noise Model for PSDF
% The TSDF needs a truncation big enough to capture the noise in the data.
% The PSDF on the other hand separates the truncation distance from the
% noise model. The probability density functions are stored as histograms
% and can therefore store complex distributions. There is no Gaussian
% assumption in the storage. This makes it possible to use a smaller
% truncation distance than the size of the noise in the data. It also means
% that there are no artifacts if a big truncation distance is used.

% Noise parameters
sigma = [0.5 3.0 0.5 2.0];
mu = [0 0 2 5];

% Metric size
trunc_dist = 2;

% PSDF quantization parameters
Np = 30;
prob_range = linspace(-trunc_dist, trunc_dist, Np);

figure(1); clf
for i=1:numel(sigma)
    subplot(2,2,i)
    plot_normpdf(mu(i), sigma(i), prob_range, [-5 5], [0 1])
end

%% Adjusting histogram interval widths to better capture wide distributions
tx = pi*0.45;
prob_range = tan(linspace(-tx, tx, Np));

% clf
% x = linspace(-pi*0.49, pi*0.49, 300);
% y = tan(x);
% plot(x, y)
% %%

figure(2); clf
for i=1:numel(sigma)
    subplot(2,2,i)
    plot_normpdf(mu(i), sigma(i), prob_range, [-5 5], [0 1])
end

%% Modeling the span of possible values
% The open and closed hypotheses give the upper and lower bound on the
% distance to the closest surface. A naive noise model is to make this
% interval fuzzy by interpreting it as a normal distribution.
%
% lo < x < hi
% lo < x & x < hi
% p(lo < x & x < hi) = p(lo < x) * p(x < hi)
% p(lo < x) = normcdf(x, lo, sigma)
% p(x < hi) = !p(hi < x) = 1 - normcdf(x, hi, sigma)
% p(lo < x & x < hi) = normcdf(x, lo, sigma) .* (1 - normcdf(x, hi, sigma))
%
% pdf = p(lo < x & x < hi) = normcdf(x, lo, sigma) .* (1 - normcdf(x, hi, sigma))
% cdf = integral(pdf, -inf, x)

lo = -2;
hi = 3;
sigma = 0.5;

% Compute histogram
prob_range_limits = prob_range(1:end-1) + diff(prob_range)*0.5;
prob_range_limits2 = [-inf, prob_range_limits, inf];

% Probability density function from the open and closed hypotheses
fun = @(x) normcdf(x, lo, sigma) .* (1 - normcdf(x, hi, sigma));

% Compute integrals over the ranges in the histogram
prob = zeros(1, Np);
for i = 1:Np
    prob(i) = integral(fun, prob_range_limits2(i), prob_range_limits2(i+1));
end

% Normalize probability
prob = prob / sum(prob);

figure(3); clf;
prob_scaled = prob./[1 diff(prob_range_limits) 1];
bar(prob_range, prob_scaled, 'FaceColor', [.5 .5 .5]);
hold on
X = linspace(-5,5,100);
plot(X, fun(X)/integral(fun, -inf, inf), 'LineWidth', 2);
xlim([-5 5])
ylim([0 0.25])
