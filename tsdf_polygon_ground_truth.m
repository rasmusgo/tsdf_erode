function [ Z, W ] = tsdf_polygon_ground_truth( N, M, trunc_dist, P, linear_weighting )
%TSDF_POLYGON Simulate 1D depth images of a polygon from arbitrary angle
%   [ Z W ] = tsdf_polygon( N, M, trunc_dist, ang, P, smoothW )
% P should be a 2 x n matrix, where n is the number of vertices 

%% Draw the polygon
% Draw one line segment at a time

img = zeros(N);
scale = 1 / (M(2)-M(1));
P2 = round((P-M(1))*scale + 1);
for i = 2:size(P2,2)
    img = draw_line(img, P2(:,i), P2(:,i-1), 1);
end

Z = bwdist(img) / scale;

% Is the polygon closed?
if P(:,1) == P(:,end)
    % Set the inside to negative distances
    mask = poly2mask(P2(1,:), P2(2,:), N, N);
    Z(mask) = -Z(mask);
    % Subtract a small offset to compensate for finite grid resolution
    Z = Z - 0.5 / scale;
end

% Compute weights from distances
W = double(Z >= -trunc_dist);
if linear_weighting
    Idx = (Z >= -trunc_dist) & (Z < 0);
    W2 = (trunc_dist + Z) / trunc_dist;
    W(Idx) = W2(Idx);
end

% Truncate the distances
Z(Z>trunc_dist) = trunc_dist;
Z(Z<-trunc_dist) = 0;

end
