%% Create 2D TSDF
N = 100;
M = linspace(-10, 10, N);
Np = 10;
r = 5;
trunc_dist = 2.5;
prob_range = linspace(-trunc_dist, trunc_dist, Np);
smoothW = true;
ang360 = [0:199] / 200 * pi*2;

%% Create ground truth of sphere (circle)
[X, Y] = meshgrid(M, M);
Z = sqrt(X.^2 + Y.^2) - r;
Z(Z>trunc_dist) = trunc_dist;
Z(Z<-trunc_dist) = -trunc_dist;

ang = linspace(0,2*pi,361);
C = [cos(ang); sin(ang)] * r;

figure(1)
plot_tsdf(Z, ones(N), trunc_dist, C)

%% Simulate 1D depth images of a circle
ang = 0/180*pi;
[Z, W] = tsdf_circle( N, M, trunc_dist, ang, 0, 0, r, smoothW );

figure(2)
plot_tsdf(Z, W, trunc_dist, C)

%% Simulate scanning a circle, two snapshots 90 degrees apart
tsdf_values = Z .* W + Z' .* W';
tsdf_weights = W + W';
tsdf_values = tsdf_values ./ tsdf_weights;
tsdf_values(tsdf_weights==0) = 0;

figure(3)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist, C)

%% Simulate scanning a circle, 90 degrees smooth turn
ang = linspace(0, pi*0.5, 100);
[tsdf_values, tsdf_weights] = tsdf_circle(N, M, trunc_dist, ang, 0, 0, r, smoothW);

figure(4)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist, C)

%% Simulate 1D depth images of a plane from arbitrary angle
P = [-r,0; r,0]';
ang = pi*-0.2;
[Z, W] = tsdf_polygon(N, M, trunc_dist, ang, P, smoothW);
plane_geometry = [-r, -1; -r, 1; -r,0; r,0; r,1; r,-1]';

figure(5)
plot_tsdf(Z, W, trunc_dist, plane_geometry)

%% Simulate scanning a plane with two snapshots
tsdf_values = Z .* W + fliplr(Z .* W);
tsdf_weights = W + fliplr(W);
tsdf_values = tsdf_values ./ tsdf_weights;
tsdf_values(tsdf_weights==0) = 0;

figure(6)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist, plane_geometry)

%% Simulate scanning a plane, turning < 90 degrees
ang = linspace(0, pi*0.45, 100);
[tsdf_values, tsdf_weights] = tsdf_polygon(N, M, trunc_dist, ang, P, smoothW);

figure(7)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist, plane_geometry)

%% Groundtruth for plane
[tsdf_values, tsdf_weights] = tsdf_polygon_ground_truth(N, M, trunc_dist, P, smoothW);

figure(8)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist, P)

%% Simulate scanning a triangle
P = [-5,-2.5; 5,-2.5; -5,2.5; -5,-2.5]';
[tsdf_values, tsdf_weights] = tsdf_polygon(N, M, trunc_dist, ang360, P, smoothW);

figure(9)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist, P)

%% Simulate scanning three triangles
P = [-4, 4; 4, 4; -4, 0; 0, 0; -4, -2; -2, -2; -4, -3; -4, 4 ]';
[tsdf_values, tsdf_weights] = tsdf_polygon(N, M, trunc_dist, ang360, P, smoothW);

figure(10)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist, P)

%% Simulate scanning L-shape
P = [-4, -4; 4, -4; 4, 0; 0, 0; 0, 4; -4, 4; -4, -4]'; % L-shape
ang = 255/180*pi;
[tsdf_values, tsdf_weights] = tsdf_polygon(N, M, trunc_dist, ang360, P, smoothW);

figure(11)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist, P)

%% Simulate scanning a cube
P = [-4, 4; 4, 4; 4, -4; -4, -4; -4, 4]'; % cube
[tsdf_values, tsdf_weights] = tsdf_polygon(N, M, trunc_dist, ang360, P, smoothW);

figure(12)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist, P)

%% Groundtruth for polygon
[tsdf_values, tsdf_weights] = tsdf_polygon_ground_truth(N, M, trunc_dist, P, smoothW);

figure(13)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist, P)

%% Simulate scanning a 'C'
P = [4.0, 0.0, -4.0, -4.0,  0.0,  4.0,  5.2,  5.5,  4.7,  3.0,  0.9, -0.1, -0.1, 0.9, 3.0, 4.7, 5.5, 5.2, 4.0; ...
     8.0, 8.0,  4.0, -4.0, -8.0, -8.0, -7.1, -4.6, -3.5, -4.0, -3.0, -1.0,  1.0, 3.0, 4.0, 3.5, 4.6, 7.1, 8.0];

[tsdf_values, tsdf_weights] = tsdf_polygon(N, M, trunc_dist, ang360, P, smoothW);

figure(14)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist, P)

%% Simulate scanning the letter 'E'
P = [-1, 0; -1, 7; 5, 7; 5, 4; 4, 5; 1, 5; 1, 1; 2, 1; 3, 2; 3, -2; 2, -1; 1, -1; 1, -5; 4, -5; 5, -4; 5, -7; -1, -7; -1, 0]';
[tsdf_values, tsdf_weights] = tsdf_polygon(N, M, trunc_dist, ang360, P, smoothW);

figure(15)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist, P)

%% Groundtruth for polygon
[tsdf_values, tsdf_weights] = tsdf_polygon_ground_truth(N, M, trunc_dist, P, smoothW);

figure(16)
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist, P)

%% Show visible lines
ang = 2/5*-90/180*pi;
L = visible_lines(ang, P);

figure(17)
plot([L(1,:); L(3,:)], [L(2,:); L(4,:)])
axis equal
axis([-10,10, -10,10])

%% Probability stuff
ang = linspace(-90, 0, 5)/180*pi;
ang = ang360;
%ang = 2/5*-90/180*pi;
%ang = -90*4/4 / 180*pi;
prob = psdf_polygon(N, M, prob_range, ang, P);

figure(18)
Pr = reshape(logspace(5, 0, Np), 1,1,Np);
prior = repmat(Pr,N);
plot_psdf(prob .* prior, prob_range, P)

figure(19)
[tsdf_values, tsdf_weights] = tsdf_polygon(N, M, trunc_dist, ang, P, smoothW);
plot_tsdf(tsdf_values, tsdf_weights, trunc_dist, P)

%% Ground truth polygon
prob = psdf_polygon_ground_truth(N, M, prob_range, P);

figure(20)
Pr = reshape(logspace(5, 0, Np), 1,1,Np);
prior = repmat(Pr,N);
plot_psdf(prob .* prior, prob_range, P)
