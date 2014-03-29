function [ W ] = erode_weights( Z, W, trunc_dist )
%ERODE_WEIGHTS Compute new weights with decreased weights near edges

%%
% Measure the difference in distance in the neighborhood of each voxel
% and use this distance as a measure of "edgyness"

%E = Z(2:end,2:end) - Z(1:end-1,2:end) - Z(2:end,1:end-1) + Z(1:end-1,1:end-1);
Ex = Z(2:end-1, 3:end) - Z(2:end-1, 1:end-2);
Ey = Z(3:end, 2:end-1) - Z(1:end-2, 2:end-1);
E = sqrt(Ex.^2 + Ey.^2);

k = fspecial('sobel');
Ex = imfilter(Z, k, 'replicate');
Ey = imfilter(Z, k', 'replicate');
E = sqrt(Ex.^2 + Ey.^2);

se = strel('square',5);
E2 = imdilate(E,se);
W2 = (20 - E2) / 10;
W2(W2 > 1) = 1;
W2(W2 < 0.1) = 0.1;
W = W2 .* W;

% clf
% subplot(2,2,1)
% imagesc( E )
% colorbar
% axis equal xy tight
% subplot(2,2,2)
% imagesc( W2 )
% colorbar
% axis equal xy tight
% subplot(2,2,3)
% imagesc( Z )
% colorbar
% axis equal xy tight
% subplot(2,2,4)
% imagesc( W )
% colorbar
% axis equal xy tight
end
