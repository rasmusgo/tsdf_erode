function [ Z W ] = tsdf_plane( N, M, trunc_dist, ang, w, smoothW )
%TSDF_PLANE Simulate 1D depth images of a plane from arbitrary angle
[X Y] = meshgrid(M, M);
k = tan(ang); m = 0;
R = [cos(ang) -sin(ang);sin(ang) cos(ang)];
X2 = R(1,1) * X + R(1,2) * Y;
Y2 = R(2,1) * X + R(2,2) * Y;

x = M;
y = k*x + m;
y(abs(x) > w*cos(ang)) = 100;
Z = zeros(N,N) - 100;
for ix = 1:N
    for iy = 1:N
        ii = round(50+5*X2(iy,ix));
        if ii > 0 && ii <= N
            Z(iy,ix) = y(ii)-Y2(iy,ix);
        end
    end
end
W = double(Z >= -trunc_dist);
if smoothW
    Idx = (Z >= -trunc_dist) & (Z < 0);
    W2 = (trunc_dist + Z) / trunc_dist;
    W(Idx) = W2(Idx);
end

Z(Z>trunc_dist) = trunc_dist;
Z(Z<-trunc_dist) = 0;
end

