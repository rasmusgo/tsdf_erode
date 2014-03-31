function [ Z W ] = tsdf_circle( N, M, trunc_dist, ang, x, y, r, smoothW )
%TSDF_CIRCLE Simulate 1D depth images of a circle from arbitrary angle

%% Call self repeatedly if multiple angles are given
if numel(ang) > 1
    Z = zeros(N);
    W = zeros(N);
    for a = ang
        [Za, Wa] = tsdf_circle(N, M, trunc_dist, a, x, y, r, smoothW);
        Z = Z + Za .* Wa;
        W = W + Wa;
    end

    Z = Z ./ W;
    Z(W==0) = 0;
    return
end

%% Create coordinate system with offset
[X, Y] = meshgrid(M-x, M-y);
R = [cos(ang) -sin(ang);sin(ang) cos(ang)];
X2 = R(1,1) * X + R(1,2) * Y;
Y2 = R(2,1) * X + R(2,2) * Y;

%% Raytrace in eye-space
% Compute projective distance
Z = -sqrt(r^2 - X2.^2) - Y2;
Z(imag(Z)~=0) = 100;
Z = real(Z);

% Compute weights
if smoothW
    W = 1 + Z / trunc_dist;
    W(Z >= 0) = 1;
    W(Z <= -trunc_dist) = 0;
else
    W = double(Z >= -trunc_dist);
end

% Truncate distances
Z(Z > trunc_dist) = trunc_dist;
Z(Z < -trunc_dist) = 0;
end

