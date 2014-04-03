function prob = psdf_polygon_ground_truth( N, M, trunc_dist, P )
%% Draw the polygon
% Draw one line segment at a time
img = zeros(N);
scale = 1 / (M(2)-M(1));
P2 = round((P-M(1))*scale + 1);
for i = 2:size(P2,2)
    img = draw_line(img, P2(:,i), P2(:,i-1), 1);
end

% Perform distance transform
Zopen = bwdist(img) / scale;

% Is the polygon closed?
if P(:,1) == P(:,end)
    % Set the inside to negative distances by swapping sign
    mask = poly2mask(P2(1,:), P2(2,:), N, N);
    Zclosed = Zopen;
    Zclosed(mask) = -Zclosed(mask);
else
    Zclosed = Zopen;
end

% Truncate the distances
Zopen(Zopen > trunc_dist) = trunc_dist;
Zopen(Zopen < -trunc_dist) = -trunc_dist;
Zclosed(Zclosed > trunc_dist) = trunc_dist;
Zclosed(Zclosed < -trunc_dist) = -trunc_dist;

%% Convert range to probability distribution
Np = 10;
%prob_range = linspace(-trunc_dist, trunc_dist, Np);
Zopen_ind = 1 + round((Np-1) * (trunc_dist + Zopen) / (trunc_dist * 2));
Zclosed_ind = 1 + round((Np-1) * (trunc_dist + Zclosed) / (trunc_dist * 2));

prob = zeros(N,N,Np);
for i = 1:Np
    prob(:,:,i) = double(i <= Zopen_ind & i >=  Zclosed_ind);
end

%plot_probability(prob)
end
