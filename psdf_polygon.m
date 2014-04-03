function prob = psdf_polygon( N, M, trunc_dist, ang, P )
%% TODO: Call self repeatedly if multiple angles are given
if numel(ang) > 1
    prob = psdf_polygon( N, M, trunc_dist, ang(1), P );
    for a = ang(2:end)
        % call self
        p = psdf_polygon( N, M, trunc_dist, a, P );
        prob = prob .* p;
    end
    return
end

%% Create rotation matrix
R = [cos(ang) -sin(ang);sin(ang) cos(ang)];

%% Find line segments visible from current view
L = visible_lines(ang, P);

%% Open hypothesis
% Draw the outline of the polygon and make a distance transform
img = zeros(N);
scale = 1 / (M(2)-M(1));
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

% Connect the lines and extend away from camera
p1 = L2(1:2, 1)   + round(1.5*N * R(2,:)');
p2 = L2(3:4, end) + round(1.5*N * R(2,:)');
shape = [p1(:); L2(:); p2]; %L2(1:2, 1), p1, p2, L2(3:4, end)];
P3 = reshape(shape, 2, numel(shape)/2);
m = [true any(P3(:,1:end-1) ~= P3(:,2:end), 1)];
P3 = P3(:,m);

% plot(P3(1,:), P3(2,:), '-o')
% axis equal
% axis([0,100, 0,100])
%%
% Draw lines with ShapeInserter
lineInserter = vision.ShapeInserter(...
    'Shape', 'Polygons', ...
    'BorderColor', 'White');
img = step(lineInserter, img, P3(:));

% Perform distance transform
Zclosed = bwdist(img) / scale;

% Set the inside to negative distances
% fillInserter = vision.ShapeInserter(...
%     'Shape', 'Polygons', ...
%     'Fill', true, ...
%     'FillColor', 'Custom', ...
%     'CustomFillColor', -1, ...
%     'Opacity', 1.0);

% Swap sign of the area behind the visible lines
%signs = step(fillInserter, ones(N)*0, P3(:));
%signs = step(lineInserter, signs, P3(:));
% imagesc(signs, [-trunc_dist, trunc_dist])
% axis xy equal tight

%Zclosed = Zclosed .* signs - 0.5 / scale;

% Set the inside to negative distances
mask = poly2mask(P3(1,:), P3(2,:), N, N);
Zclosed(mask) = -Zclosed(mask);

% imagesc(Zclosed, [-trunc_dist, trunc_dist])
% axis xy equal tight

%% Truncate the distances
Zopen(Zopen > trunc_dist) = trunc_dist;
Zopen(Zopen < -trunc_dist) = -trunc_dist;
Zclosed(Zclosed > trunc_dist) = trunc_dist;
Zclosed(Zclosed < -trunc_dist) = -trunc_dist;

%% Convert range to probability distribution
Np = 10;
% Zopen_ind = 1 + round((trunc_dist + Zopen) * (Np-1) / (trunc_dist * 2));
% Zclosed_ind = 1 + round((trunc_dist + Zclosed) * (Np-1) / (trunc_dist * 2));

prob_range = linspace(-trunc_dist, trunc_dist, Np);
epsilon = 2/scale;
prob = zeros(N,N,Np);
for i = 1:Np
    z = prob_range(i);
    prob(:,:,i) = double(z <= (Zopen+epsilon) & z >=  (Zclosed-epsilon));
end

%plot_probability(prob)
end
