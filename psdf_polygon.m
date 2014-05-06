function prob = psdf_polygon( N, M, prob_range, ang, P, multiplier )
if nargin < 6
    multiplier = 1;
end
%% Call self repeatedly if multiple angles are given
if numel(ang) > 1
    prob = psdf_polygon( N, M, prob_range, ang(1), P, multiplier );
    for a = ang(2:end)
        % call self
        p = psdf_polygon( N, M, prob_range, a, P, multiplier );
        prob = prob .* p;
    end
    return
end

%% Find line segments visible from current view
L = visible_lines(ang, P);

%% Open hypothesis
% Draw the outline of the polygon and make a distance transform
%multiplier = 4;
N2 = N*multiplier;
img = zeros(N2);
scale = multiplier / (M(2)-M(1));
L2 = round((L-M(1))*scale + 1);
for i = 1:size(L2,2)
    img = draw_line(img, L2(1:2,i), L2(3:4,i), 1);
end

% Perform distance transform
Zopen = bwdist(img) / scale;

% imagesc(Zopen, [-trunc_dist, trunc_dist])
% axis xy equal tight

%% Closed hypothesis
% Draw the outline of the polygon, make a distance transform
% Fill the polygon, use this as a mask to invert the distances

% Connect the lines and extend edges away from camera
R = [cos(ang) -sin(ang);sin(ang) cos(ang)];
p1 = L2(1:2, 1)   + round(1.5*N2 * R(2,:)');
p2 = L2(3:4, end) + round(1.5*N2 * R(2,:)');
shape = [p1(:); L2(:); p2];
P3 = reshape(shape, 2, numel(shape)/2);

% Remove doubles created from already connected lines
P3 = P3(:, [true any(P3(:,1:end-1) ~= P3(:,2:end), 1)]);

% plot(P3(1,:), P3(2,:), '-o')
% axis equal
% axis([0,100, 0,100])

%% Create mask of perimeter
img = draw_polygon(img, P3, 1);

% Perform distance transform
Zclosed = bwdist(img) / scale;

% Set the inside to negative distances
mask = poly2mask(P3(1,:), P3(2,:), N2, N2);
Zclosed(mask) = -Zclosed(mask);

% imagesc(Zclosed, [prob_range(1), prob_range(end)])
% axis xy equal tight

% Undo scale multiplier
Zopen = Zopen(1:multiplier:end, 1:multiplier:end);
Zclosed = Zclosed(1:multiplier:end, 1:multiplier:end);

%% Convert range of possible values to probability distribution
epsilon = 0.5 * (M(2)-M(1));
prob_range_mid = prob_range(1:end-1) + diff(prob_range)*0.5;
prob_range_lo = [-inf, prob_range_mid] - epsilon;
prob_range_hi = [prob_range_mid, inf] + epsilon;
Np = numel(prob_range);
prob = zeros(N,N,Np);
for i = 1:Np
    prob(:,:,i) = double(Zopen >= prob_range_lo(i) & Zclosed <= prob_range_hi(i));
end

%plot_probability(prob)
end
