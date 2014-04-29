%%
clf
X = -10:0.01:10; %linspace(-10, 10, 1000);

% observation parameters
x = 1;
lo1 = -4;
hi1 = 4;
lo2 = 1;
hi2 = 1;
sigma = 1.5;
delta = 0;%1*sigma;

% Probability density function from the open and closed hypotheses
fun1 = @(x) normcdf(x, lo1-delta, sigma);
fun2 = @(x) (1 - normcdf(x, hi1+delta, sigma));
fun3 = @(x) normcdf(x, lo2-delta, sigma);
fun4 = @(x) (1 - normcdf(x, hi2+delta, sigma));

%Y = normpdf(X, 0, 1);
P_outlier = 1e-2;
W_open = 0.25;
W_closed = 0.75;
W_between = 0.75;
Z1 = W_between * fun1(X) .* fun2(X) + W_closed * exp(-0.5 * ((X - lo1)./sigma).^2) + W_open * exp(-0.5 * ((X - hi1)./sigma).^2) + P_outlier;
Z2 = W_between * fun3(X) .* fun4(X) + W_closed * exp(-0.5 * ((X - lo2)./sigma).^2) + W_open * exp(-0.5 * ((X - hi2)./sigma).^2) + P_outlier;
Z1 = Z1 / max(Z1);
Z2 = Z2 / max(Z2);
Y = Z1 .^ 1 .* Z2 .^ 1;
Y = Y / max(Y);

hold on
plot(X, [Z1; Z2; Y], 'LineWidth', 2)
plot(X, [-log(Z1); -log(Z2); -log(Y)], ':', 'LineWidth', 1)

E_Z1 = X * Z1' / sum(Z1);
E_Z2 = X * Z2' / sum(Z2);
E_Y = X * Y' / sum(Y);

[Mv_Z1, Mi_Z1] = max(Z1);
[Mv_Z2, Mi_Z2] = max(Z2);
[Mv_Y, Mi_Y] = max(Y);
M_Z1 = X(Mi_Z1);
M_Z2 = X(Mi_Z2);
M_Y = X(Mi_Y);

plot([M_Z1 M_Z1; M_Z2 M_Z2; M_Y M_Y]', [-1 2; -1 2; -1 2]', '--', 'LineWidth', 1)

legend 'Z1' 'Z2' 'Y' '-log(Z1)' '-log(Z2)' '-log(Y)' 'MAP(Z1)' 'MAP(Z2)' 'MAP(Y)'
ylim([-0.1,1.5])
