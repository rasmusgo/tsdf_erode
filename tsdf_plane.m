function [ Z W ] = tsdf_plane( N, M, trunc_dist, ang, x, y, w, smoothW )
%TSDF_PLANE Simulate 1D depth images of a plane from arbitrary angle

%% Create rotated coordinate system
[X, Y] = meshgrid(M-x, M-y);
R = [cos(ang) -sin(ang);sin(ang) cos(ang)];
X2 = R(1,1) * X + R(1,2) * Y;
Y2 = R(2,1) * X + R(2,2) * Y;

%% Raytrace in eye-space
k = tan(ang); m = 0;
Z = k*X2 + m - Y2;
Z(abs(X2) > w*cos(ang)) = trunc_dist;

W = double(Z >= -trunc_dist);
if smoothW
    Idx = (Z >= -trunc_dist) & (Z < 0);
    W2 = (trunc_dist + Z) / trunc_dist;
    W(Idx) = W2(Idx);
end

Z(Z>trunc_dist) = trunc_dist;
Z(Z<-trunc_dist) = 0;
end
