function logodds = occu_rayprob( N, M, ang, points, sigma, P_outlier )

% TODO: Project each voxel onto depthmap and use that distance to compute
% inverse measurement model.

R = [cos(ang) -sin(ang);sin(ang) cos(ang)];
points_cam = R * points;

% Start with complete uncertainty
logodds = log(ones(N));

% Compute position for each voxel in camera space
[X, Y] = meshgrid(M, M);
XY_cam = R * [X(:) Y(:)]';

% Find linear mapping from cam to depth map (points)
x = points_cam(1,:);
dx = diff(x);
dx = mean(dx(dx==dx));
I = 1:size(x, 2);
x0 = mean(x(x==x) - I(x==x)*dx);

% Find associated point in depth map for each voxel
Ix = round((XY_cam(1,:) - x0) ./ dx);
nXY = size(XY_cam, 2);
dY = nan(1, nXY);
valid_ind = Ix >= 1 & Ix <= size(points_cam, 2);
dY(valid_ind) = XY_cam(2, valid_ind) - points_cam(2, Ix(valid_ind));

% Sensor model, dY is distance to surface, negative is closer to viewer
logodds(dY < 0) = log(0.1);
logodds(0 < dY & dY < 0.5) = log(9);

end
