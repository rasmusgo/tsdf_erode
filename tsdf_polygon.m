function [ Z, W ] = tsdf_polygon( N, M, trunc_dist, ang, P, linear_weighting )
%TSDF_POLYGON Simulate 1D depth images of a polygon from arbitrary angle
%   [ Z W ] = tsdf_polygon( N, M, trunc_dist, ang, P, smoothW )
% P should be a 2 x n matrix, where n is the number of vertices 

%% Call self repeatedly if multiple angles are given
if numel(ang) > 1
    Z = zeros(N);
    W = zeros(N);
    for a = ang
        [Za, Wa] = tsdf_polygon(N, M, trunc_dist, a, P, linear_weighting);
        Z = Z + Za .* Wa;
        W = W + Wa;
    end

    Z = Z ./ W;
    Z(W==0) = 0;
    return
end

%% Create rotated coordinate system
[X, Y] = meshgrid(M, M);
R = [cos(ang) -sin(ang);sin(ang) cos(ang)];
X2 = R(1,1) * X + R(1,2) * Y;
Y2 = R(2,1) * X + R(2,2) * Y;

% Rotate the polygon
P2 = R * P;

%% Raytrace in eye-space
% Try one line segment at a time
% Keep the closest (lowest value)

Z = ones(N)*trunc_dist;
for i = 2:size(P2,2)
    dp = P2(:,i) - P2(:,i-1);
    dp = dp / norm(dp);
    k = dp(2)/dp(1);
    m = P2(2,i) - k*P2(1,i);
    lo = min( P2(1,i), P2(1,i-1) );
    hi = max( P2(1,i), P2(1,i-1) );
    Zi = k*X2 + m - Y2;
    Zi(X2 < lo | X2 > hi) = trunc_dist;
    Z = min(Z, Zi);
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
