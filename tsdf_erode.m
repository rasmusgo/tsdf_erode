%% Create 2D TSDF
width = 100;
height = 100;
tsdf_values = zeros(width, height);
tsdf_weights = zeros(width, height);

%% Create ground truth of sphere (circle)
[X Y] = meshgrid(linspace(-10, 10, 100), linspace(-10, 10, 100));
Z = sqrt(X.^2 + Y.^2) - 5;
axis equal
imagesc(Z)
colorbar
hold on
contour(Z, [0 0], '-k')

%% Simulate 1D depth images of a sphere (circle)
depth = zeros(height,1);
X = linspace(-10, 10, 100);

%% Integrate TSDF

%% Show zero crossings
contour(tsdf_values, [0 0])
