%% PSDF Tutorial
% The Probabalistic Signed Distance Function is similar to a TSDF but
% instead of having a distance and a weight for each voxel, they each
% contain a probability distribution of the possible values for the
% distance to the closest surface.
%
% Using a probability distribution makes uncertain areas explicitly
% represented. The area behind the observed geometry gets a wide
% distribution instead of a single value.
% 
% The probability distribution is updated for each frame. Two distance
% transforms are used to calculate the upper and lower limit for the
% distance to the nearest surface based on the observed surfaces. The 
% interval between the upper and lower limit is interpreted as a uniform
% probability distribution. The final PSDF volume is the product of
% all the frames probability distributions.
% 

%% Global parameters

% Volume resolution
N = 100;

% Metric size
M = linspace(-10, 10, N);
trunc_dist = 2.5;

% PSDF quantization parameters
Np = 10;
prob_range = linspace(-trunc_dist, trunc_dist, Np);
linear_weighting = true;

% Observation angles
Nang360 = 200;
ang360 = (0:(Nang360-1)) / Nang360 * pi*2;

% Radius of sphere
r = 5;

%% Input observation
% Looking at the letter 'E' from the bottom right with an orthographic
% camera gives an incomplete view of the geometry.

% The letter 'E' as a polygon
P = [-1, 0; -1, 7; 5, 7; 5, 4; 4, 5; 1, 5; 1, 1; 2, 1; 3, 2; 3, -2; 2, -1; 1, -1; 1, -5; 4, -5; 5, -4; 5, -7; -1, -7; -1, 0]';

% Compute visible line segments
ang = -45/180*pi;
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
% The closed hypothesis is formed by connecting the observed geomtry and
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
img = insertShape(img, 'Polygon', P3(:)', 'SmoothEdges', false, 'Color', 'White', 'Opacity', 1.0);
img = img(:,:,1);

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
epsilon = 0.5/scale;
prob_range_mid = prob_range(1:end-1) + diff(prob_range)*0.5;
prob_range_lo = [-inf, prob_range_mid] - epsilon;
prob_range_hi = [prob_range_mid, inf] + epsilon;
Np = numel(prob_range);
prob = zeros(N,N,Np);
for i = 1:Np
    prob(:,:,i) = double(Zopen > prob_range_lo(i) & Zclosed < prob_range_hi(i));
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

% Compute expected values
Z = sum(bsxfun(@times, reshape(prob_range, 1,1,Np), prob), 3) ./ sum(prob, 3);

figure(3); clf; hold on
imagesc(Z, [prob_range(1) prob_range(end)])
colorbar
colormap hot

% Plot ground truth polygon and observed line segments
plot(P(1,:)*5+50.5, P(2,:)*5+50.5, 'k:', 'LineWidth', 1)
plot([L(1,:); L(3,:)]*5+50.5, [L(2,:); L(4,:)]*5+50.5, 'k', 'LineWidth', 3)

axis equal xy tight

title 'Expected distance to closest surface'

samples = [30,70; 30,30; 70,50]';
short_descriptions = {'1)', '2)', '3)'};
descriptions = {'1) Uncertain', '2) Limited range', '3) Certain'};

% Mark samples in figure
for i = 1:size(samples,2)
    text(samples(1,i), samples(2,i), short_descriptions(i))
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
plot_psdf(psdf1, prob_range, P)
title 'View 1'

subplot(1,2,2)
plot_psdf(psdf2, prob_range, P)
title 'View 2'

figure(6); clf
plot_psdf(psdf3, prob_range, P)
title 'Combined views'

%%
% A result close to ground truth can be obtained by combining lots of views
% and a prior that prefers negative distances.
prob = psdf_polygon(N, M, prob_range, ang360, P);
Pr = reshape(logspace(0, -5, Np), 1,1,Np);
prior = repmat(Pr,N);

figure(7); clf
subplot(1,2,1)
plot_psdf(prob, prob_range, P)
title 'Full circle scan'

subplot(1,2,2)
plot_pdf(Pr, prob_range)
title 'Prior distribution'

figure(8); clf
plot_psdf(prob .* prior, prob_range, P)
title 'Prior applied to scan'
