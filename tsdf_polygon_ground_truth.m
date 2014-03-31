function [ Z W ] = tsdf_polygon_ground_truth( N, M, trunc_dist, P, smoothW )
%TSDF_POLYGON Simulate 1D depth images of a polygon from arbitrary angle
%   [ Z W ] = tsdf_polygon( N, M, trunc_dist, ang, P, smoothW )
% P should be a 2 x n matrix, where n is the number of vertices 

%% Create coordinate system
[X, Y] = meshgrid(M, M);

%% Draw the polygon
% Draw one line segment at a time

img = zeros(N);
scale = 1 / (M(2)-M(1));
P2 = round((P-M(1))*scale + 1);
for i = 2:size(P2,2)
    img = draw_line(img, P2(:,i), P2(:,i-1), 1);
end

Z = bwdist(img) / scale;

% Set the inside to negative distances
shapeInserter = vision.ShapeInserter(...
    'Shape', 'Polygons', ...
    'Fill', true, ...
    'FillColor', 'Custom', ...
    'CustomFillColor', -1, ...
    'Opacity', 1.0);

signs = step(shapeInserter, ones(N), P2(:));
Z = Z .* signs - 0.5 / scale;

% Compute weights from distances
W = double(Z >= -trunc_dist);
if smoothW
    Idx = (Z >= -trunc_dist) & (Z < 0);
    W2 = (trunc_dist + Z) / trunc_dist;
    W(Idx) = W2(Idx);
end

% Truncate the distances
Z(Z>trunc_dist) = trunc_dist;
Z(Z<-trunc_dist) = 0;

end