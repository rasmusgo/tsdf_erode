%% Create 2D TSDF
N = 100;
M = linspace(-10, 10, N);
r = 5;
trunc_dist = 4;
smoothW = true;
tsdf_values = zeros(N);
tsdf_weights = zeros(N);

%% Create ground truth of sphere (circle)
[X Y] = meshgrid(M, M);
Z = sqrt(X.^2 + Y.^2) - r;
Z(Z>trunc_dist) = trunc_dist;
Z(Z<-trunc_dist) = -trunc_dist;

figure(1)
plot_tsdf(Z, ones(N), trunc_dist)

%% Simulate 1D depth images of a sphere (circle)
ang = 0;
[Z, W] = tsdf_circle( N, M, trunc_dist, ang, r, smoothW );

figure(2)
plot_tsdf(Z, W, trunc_dist)

%% Integrate TSDF
tsdf_values = Z .* W + Z' .* W';
tsdf_weights = W + W';
tsdf_values = tsdf_values ./ tsdf_weights;
tsdf_values(tsdf_weights==0) = 0;

figure(3)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist)

%% Simulate 1D depth images of a plane
ang = pi*-0.25; k = tan(ang); m = 0;
depth = zeros(N,1);
x = M;
y = k*x + m;
y(abs(x) > r*cos(ang)) = 100;
[X, Y] = meshgrid(M, M);
Z = bsxfun(@minus, y, Y);
W = Z >= -trunc_dist;
Z(Z>trunc_dist) = trunc_dist;
Z(Z<-trunc_dist) = 0;

figure(4)
plot_tsdf(Z, W, trunc_dist)

%% Integrate TSDF
tsdf_values = Z .* W + Z' .* W';
tsdf_weights = W + W';
tsdf_values = tsdf_values ./ tsdf_weights;
tsdf_values(tsdf_weights==0) = 0;

figure(5)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist)

%% Simulate 1D depth images of a plane from arbitrary angle
[X Y] = meshgrid(M, M);
ang = pi*-0.2; k = tan(ang); m = 0;
[Z, W] = tsdf_plane(N, M, trunc_dist, ang, r, smoothW);

figure(6)
plot_tsdf(Z, W, trunc_dist)

%% Call tsdf_plane repeatedly and accumulate
tsdf_values = zeros(N);
tsdf_weights = zeros(N);
for ang = linspace(0, pi*0.45, 100)
    [Z, W] = tsdf_plane(N, M, trunc_dist, ang, r, smoothW);
    tsdf_values = tsdf_values + Z .* W;
    tsdf_weights = tsdf_weights + W;
end

tsdf_values = tsdf_values ./ tsdf_weights;
tsdf_values(tsdf_weights==0) = 0;

figure(7)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist)
hold on
plot([r r]*5+50, [-r r]*5+50, '-k')
