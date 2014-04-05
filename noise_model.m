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

