function prob = psdf_polygon_ground_truth( N, M, prob_range, P )
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

%% Convert range of possible values to probability distribution
epsilon = 0.5/scale;
prob_range_mid = prob_range(1:end-1) + diff(prob_range)*0.5;
prob_range_lo = [-inf, prob_range_mid] - epsilon;
prob_range_hi = [prob_range_mid, inf] + epsilon;
Np = numel(prob_range);
prob = zeros(N,N,Np);
for i = 1:Np
    prob(:,:,i) = double(Zopen > prob_range_lo(i) & Zclosed < prob_range_hi(i));
end

%plot_probability(prob)
end
