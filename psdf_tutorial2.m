%% PSDF Tutorial 2
% Replacing histogram ranges with samplings of a pdf

%% Global parameters

% Volume resolution
N = 100;

% Metric size
M = linspace(-10, 10, N);
trunc_dist = 2.0;

% PSDF quantization parameters
Np = 20;
prob_range = linspace(-trunc_dist, trunc_dist, Np);

% Observation angles
Nang360 = 200;
ang360 = (0:(Nang360-1)) / Nang360 * pi*2;

% Observation uncertainty
sigma_x = 0.1;
sigma_y = 0.1;
sigma = 0.25;

%% Input observation
% Looking at the letter 'E' from the bottom right with an orthographic
% camera gives an incomplete view of the geometry.

% The letter 'E' as a polygon
P = [-1, 0; -1, 7; 5, 7; 5, 4; 4, 5; 1, 5; 1, 1; 2, 1; 3, 2; 3, -2; 2, -1; 1, -1; 1, -5; 4, -5; 5, -4; 5, -7; -1, -7; -1, 0]';

P = [cos(linspace(-pi,pi, 10)) ; sin(linspace(-pi,pi, 10))]*5;

% Compute visible line segments
ang = -30/180*pi;
L = visible_lines(ang, P);

figure(1)
clf; hold on
title 'Line segments visible from bottom right'
plot(P(1,:), P(2,:), 'k:', 'LineWidth', 1)
plot([L(1,:); L(3,:)], [L(2,:); L(4,:)], 'k', 'LineWidth', 3)
axis('square', [-10,10, -10,10])

%% Open and closed hypotheses
% The area behind the visible line segments is unknown. Any geometry could
% exist there. The two extremes are 1) Nothing exists except what we
% observe, 2) Everything behind the observed lines is solid. The first case
% is called the _open hypothesis_ and the second case is called the _closed
% hypothesis_. These two hypotheses together gives the range of possible
% distances to the closest surface for each voxel.

%%
% The open hypothesis is formed by making a distance transform of the
% observed geometry.

% Create mask of observed geometry
img = zeros(N);
scale = 1 / (M(2)-M(1));
L2 = round((L-M(1))*scale + 1);
for i = 1:size(L2,2)
    img = draw_line(img, L2(1:2,i), L2(3:4,i), 1);
end

% Perform distance transform
Zopen = bwdist(img) / scale;

%%
% The closed hypothesis is formed by connecting the observed geometry and
% extruding it backwards, creating a perimeter around the area that could
% potentially be solid. A distance transform is performed on this
% perimeter. The values on the outside are correct but the values on the
% inside must be negated to label them as being inside a solid area.

% Connect the lines and extend edges away from camera
R = [cos(ang) -sin(ang);sin(ang) cos(ang)];
p1 = L2(1:2, 1)   + round(1.5*N * R(2,:)');
p2 = L2(3:4, end) + round(1.5*N * R(2,:)');
shape = [p1(:); L2(:); p2];
P3 = reshape(shape, 2, numel(shape)/2);

% Remove doubles created from already connected lines
P3 = P3(:, [true any(P3(:,1:end-1) ~= P3(:,2:end), 1)]);

% Create mask of perimeter
img = draw_polygon(img, P3, 1);

% Perform distance transform
Zclosed = bwdist(img) / scale;

% Set the inside to negative distances
mask = poly2mask(P3(1,:), P3(2,:), N, N);
Zclosed(mask) = -Zclosed(mask);

%%
% The distances are truncated and color coded to ease visualization.

figure(2)
plot_open_closed(Zopen, Zclosed, trunc_dist, P)

%% Creating the probability distributions
% The open and closed hypotheses give the maximal and minimal possible
% distance to the closest surface for each voxel. This range is stored as a
% histogram over possible distances. Each voxel has its own histogram. The
% probability density functions are not normalized and do not need to be.

% Convert range of possible values to probability distribution
epsilon = 0.5 * (M(2)-M(1));
prob_range_mid = prob_range(1:end-1) + diff(prob_range)*0.5;
prob_range_lo = [-inf, prob_range_mid] - epsilon;
prob_range_hi = [prob_range_mid, inf] + epsilon;
Np = numel(prob_range);
prob = zeros(N,N,Np);

% TODO: Replace sharp "perfect" distribution with uncertain distribution
%       to model noise in the observations
for i = 1:Np
    prob(:,:,i) = double(Zopen >= prob_range_lo(i) & Zclosed <= prob_range_hi(i));
end

%% Analyzing the probability distributions
% The expected distance to the closest surface can be misleading. The
% expected value will be zero behind much of the observed line segments
% from a single view. This is because the possible range of values in these
% areas are symmetric around zero. The area closely in front of the visible
% line segments have certain values because the open and closed hypotheses
% agrees on the distance to the closest surface. There are also areas where
% the uncertainty is limited. These areas are closer to the border of the
% observed area than to any observed surface.
%
% A better measurement is the value with maximum likelihood with a minimal
% prior towards lower values.

% Compute expected values
Z = sum(bsxfun(@times, reshape(prob_range, 1,1,Np), prob), 3) ./ sum(prob, 3);
[m, ind] = max(prob, [], 3);
Z2 = reshape(prob_range(reshape(ind, numel(ind), 1)), size(ind));

figure(3); clf; hold on
imagesc(Z2, [prob_range(1) prob_range(end)])
colorbar
colormap hot

% Plot ground truth polygon and observed line segments
plot(P(1,:)*5+50.5, P(2,:)*5+50.5, 'k:', 'LineWidth', 1)
plot([L(1,:); L(3,:)]*5+50.5, [L(2,:); L(4,:)]*5+50.5, 'k', 'LineWidth', 3)

axis equal xy tight

title 'Expected distance to closest surface'

samples = [30,70; 30,30; 70,50]';
short_descriptions = {'(1)', '(2)', '(3)'};
descriptions = {'(1) Uncertain', '(2) Limited range', '(3) Certain'};

% Mark samples in figure
for i = 1:size(samples,2)
    text(samples(1,i), samples(2,i), short_descriptions(i), 'Background', [1,1,1])
end

% Show histograms for selected samples
figure(4); clf
for i = 1:size(samples,2)
    subplot(1,size(samples,2),i)
    plot_pdf(prob(samples(2,i),samples(1,i),:), prob_range);
    title(descriptions(i))
end

%% Combining probability distributions
% Probability distributions from different views can be combined simply by
% multiplication.
psdf1 = prob;
psdf2 = psdf_polygon(N, M, prob_range, 0, P);
psdf3 = psdf1 .* psdf2;

figure(5); clf
subplot(1,2,1)
plot_psdf_max(psdf1, prob_range, P)
title 'View 1'

subplot(1,2,2)
plot_psdf_max(psdf2, prob_range, P)
title 'View 2'

figure(6); clf
plot_psdf_max(psdf3, prob_range, P)
title 'Combined views'

%%
% A result close to ground truth can be obtained by combining lots of
% views. There is no need to use a strong prior when using maximimum
% likelihood instead of expected value. 
prob = psdf_polygon(N, M, prob_range, ang360, P);

figure(7); clf
plot_psdf_max(prob, prob_range, P)
title 'Maximal likelihood values of 360 scan'

%% Probability distribution of ray with uncertainty
% A single depth measurement can be interpreted as "The line from the
% camera to the observed point contains open space and the observed point
% is on a surface". The open hypotheses is the distance to the observed
% point, all that is not observed as matter is considered to be open space.
% The closed hypotheses is the negative distance to the line, all that is
% not observed to be open space is considered to be solid matter.

% Define ray
ray = [3 0.001 10 -7]';

% Find line segments visible from current view
%L = visible_lines(ang, P);

% Create mask of observed geometry
img = zeros(N);
scale = 1 / (M(2)-M(1));
ray2 = round((ray-M(1))*scale + 1);

% Compute distance to observed point
img(ray2(2), ray2(1)) = 1;
Zopen = bwdist(img) / scale;

% Compute distance to observed free space
img = draw_line(img, ray2(1:2), ray2(3:4), 1);
Zclosed = -bwdist(img) / scale;

% Plot open and closed hypotheses
figure(8); clf
plot_open_closed(Zopen, Zclosed, trunc_dist, P);

%% Model uncertainty in the measurement
P_outlier = 1e-2;
W_open = 0.15;
W_closed = 0.16;
W_between = 0.75;
fun1 = @(x, lo, hi, sigma) W_between * normcdf(x, lo, sigma) .* (1 - normcdf(x, hi, sigma)) + W_closed * exp(-0.5 * ((x - lo)./sigma).^2) + W_open * exp(-0.5 * ((x - hi)./sigma).^2) + P_outlier;
fun2 = @(x1, x2, sigma) W_between * normcdf(x1, 0, sigma) .* (1 - normcdf(x2, 0, sigma)) + W_closed * exp(-0.5 * (x1./sigma).^2) + W_open * exp(-0.5 * (x2./sigma).^2) + P_outlier;

X = linspace(-10, 10, 100);
Y1 = fun1(X, -1, 3, sigma);
Y2 = fun2(X+1, X-3, sigma);
Y3 = fun2(reshape(X+1, 10, 10), reshape(X-3, 10, 10), sigma);
Y3 = Y3(:)';

figure(9); clf
plot(X,[Y1;Y2;Y3])
%%
% Apply uncertainty to casted ray

%dists = repmat(reshape(prob_range, 1, 1, Np), N, N);

for i = 1:Np
    prob(:,:,i) = double(fun2(prob_range(i) - Zclosed, prob_range(i) - Zopen, sigma));
end

% "Normalize" pdf
prob = bsxfun(@rdivide, prob, max(prob, [], 3));

figure(10); clf
plot_psdf_max(prob, prob_range, P)
title 'Maximal likelihood values of single point observation'

samples = [75,45; ray2(1), ray2(2); 60,50]';
short_descriptions = {'(1)', '(2)', '(3)'};
descriptions = {'(1) Free space', '(2) Observed point', '(3) Behind surface'};

% Mark samples in figure
for i = 1:size(samples,2)
    text(samples(1,i), samples(2,i), short_descriptions(i), 'Background', [1,1,1])
end

% Show histograms for selected samples
figure(11); clf
for i = 1:size(samples,2)
    subplot(1,size(samples,2),i)
    plot_pdf(prob(samples(2,i),samples(1,i),:), prob_range);
    title(descriptions(i))
end

%% Raytrace against geometry
Nrays = 100;
points = raytrace(ang, P, Nrays, sigma_x, sigma_y);

figure(12); clf; hold on;
title 'Intersection points with the polygon'
plot(P(1,:), P(2,:), 'k:', 'LineWidth', 1)
plot([L(1,:); L(3,:)], [L(2,:); L(4,:)], 'k', 'LineWidth', 3)
plot(points(1,:), points(2,:), 'xr', 'LineWidth', 2)
axis('square', [-10,10, -10,10])

%%
% Compute probability distribution for each ray separately
prob = psdf_rayprob( N, M, prob_range, ang, points, sigma, P_outlier, W_open, W_closed, W_between );

figure(13); clf
plot_psdf_max(prob, prob_range, P)
title 'Maximal likelihood values of noisy measurement'

%%
prob = psdf_raytracepolygon( N, M, prob_range, ang360(1:20:end), P, Nrays, sigma_x, sigma_y, sigma, P_outlier, W_open, W_closed, W_between );

%
figure(14); clf
plot_psdf_max(prob, prob_range, P)
title 'Maximal likelihood values of noisy measurement'

samples = [75,45; 66,50; 70,20]';
short_descriptions = {'(1)', '(2)', '(3)'};
descriptions = {'(1) Free space', '(2) Observed point', '(3) Behind surface'};

% Mark samples in figure
for i = 1:size(samples,2)
    text(samples(1,i), samples(2,i), short_descriptions(i), 'Background', [1,1,1])
end

% Show histograms for selected samples
figure(15); clf
for i = 1:size(samples,2)
    subplot(1,size(samples,2),i)
    plot_pdf(prob(samples(2,i),samples(1,i),:), prob_range);
    title(descriptions(i))
end

%% Occupancy grid
logodds = occu_raytracepolygon( N, M, ang360(1:1:end), P, Nrays, sigma_x, sigma_y, sigma, P_outlier );

%
figure(14); clf
plot_occu(logodds, P)
