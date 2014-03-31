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
ang = 0/180*pi;
[Z, W] = tsdf_circle( N, M, trunc_dist, ang, 0, 0, r, smoothW );

figure(2)
plot_tsdf(Z, W, trunc_dist)

%% Integrate TSDF
tsdf_values = Z .* W + Z' .* W';
tsdf_weights = W + W';
tsdf_values = tsdf_values ./ tsdf_weights;
tsdf_values(tsdf_weights==0) = 0;

figure(3)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist)

%% Call tsdf_circle repeatedly and accumulate
tsdf_values = zeros(N);
tsdf_weights = zeros(N);
for ang = linspace(0, pi*0.5, 100)
    [Z, W] = tsdf_circle(N, M, trunc_dist, ang, 0, 0, r, smoothW);
    tsdf_values = tsdf_values + Z .* W;
    tsdf_weights = tsdf_weights + W;
end

tsdf_values = tsdf_values ./ tsdf_weights;
tsdf_values(tsdf_weights==0) = 0;

figure(4)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist)

%% Simulate 1D depth images of a plane from arbitrary angle
[X Y] = meshgrid(M, M);
ang = pi*-0.2; k = tan(ang); m = 0;
[Z, W] = tsdf_plane(N, M, trunc_dist, ang, 0, 0, r, smoothW);

figure(5)
plot_tsdf(Z, W, trunc_dist)

%% Integrate TSDF
tsdf_values = Z .* W + fliplr(Z .* W);
tsdf_weights = W + fliplr(W);
tsdf_values = tsdf_values ./ tsdf_weights;
tsdf_values(tsdf_weights==0) = 0;

figure(6)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist)

%% Call tsdf_plane repeatedly and accumulate
tsdf_values = zeros(N);
tsdf_weights = zeros(N);
for ang = linspace(0, pi*0.45, 100)
    [Z, W] = tsdf_plane(N, M, trunc_dist, ang, 0, 0, r, smoothW);
    tsdf_values = tsdf_values + Z .* W;
    tsdf_weights = tsdf_weights + W;
end

tsdf_values = tsdf_values ./ tsdf_weights;
tsdf_values(tsdf_weights==0) = 0;

figure(7)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist)
hold on
plot([r r]*5+50.5, [-r r]*5+50.5, '-k')

%% Call tsdf_polygon repeatedly and accumulate
P = [-5,-2.5; 5,-2.5; -5,2.5; -5,-2.5]';
tsdf_values = zeros(N);
tsdf_weights = zeros(N);
for ang = linspace(0, pi*2, 100)
    [Z, W] = tsdf_polygon(N, M, trunc_dist, ang, P, smoothW);
    W = erode_weights(Z,W,trunc_dist);
    tsdf_values = tsdf_values + Z .* W;
    tsdf_weights = tsdf_weights + W;
end

tsdf_values = tsdf_values ./ tsdf_weights;
tsdf_values(tsdf_weights==0) = 0;

figure(8)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist)
hold on
plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-k', 'LineWidth', 3)
plot(P(1,:)*5+50.5, P(2,:)*5+50.5, '-w', 'LineWidth', 1)
