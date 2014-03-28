function [ Z W ] = tsdf_circle( N, M, trunc_dist, ang, r, smoothW )
%TSDF_PLANE Simulate 1D depth images of a plane from arbitrary angle
%   View angle is not implemented yet
x = M;
y = -sqrt(r^2 - x.^2);
y(imag(y)~=0) = 100;
y = real(y);
[X, Y] = meshgrid(M, M);
Z = bsxfun(@minus, y, Y);
if smoothW
    W = 1 + Z / trunc_dist;
    W(Z >= 0) = 1;
    W(Z <= -trunc_dist) = 0;
else
    W = double(Z >= -trunc_dist);
end

Z(Z > trunc_dist) = trunc_dist;
Z(Z < -trunc_dist) = 0;
end

