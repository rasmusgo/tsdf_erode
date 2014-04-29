function prob = psdf_rayprob( N, M, prob_range, ang, points, sigma, P_outlier, W_open, W_closed, W_between )

%% Call self repeatedly if multiple points are given
if size(points, 2) > 1
    prob = psdf_rayprob( N, M, prob_range, ang, points(:,1), sigma, P_outlier, W_open, W_closed, W_between );
    for i = 2:size(points, 2)
        % call self
        p = psdf_rayprob( N, M, prob_range, ang, points(:,i), sigma, P_outlier, W_open, W_closed, W_between );
        prob = prob .* p;
    end
    return
end

%% No points and NaN values give no information
if size(points, 2) == 0 || any(isnan(points))
    prob = ones(N,N,numel(prob_range));
    return;
end

%% Create ray from point
ray = [points; points + [-sin(ang); -cos(ang)] * (M(end) - M(1))*1.5];

%% Create mask of observed geometry
img = zeros(N);
scale = 1 / (M(2)-M(1));
ray2 = round((ray-M(1))*scale + 1);

% Compute distance to observed point
img(ray2(2), ray2(1)) = 1;
Zopen = bwdist(img) / scale;

% Compute distance to observed free space
img = draw_line(img, ray2(1:2), ray2(3:4), 1);
Zclosed = -bwdist(img) / scale;

%% Model uncertainty in the measurement
%fun2 = @(x1, x2, sigma) normcdf(x1, 0, sigma) .* (1 - normcdf(x2, 0, sigma));
fun2 = @(x1, x2, sigma) W_between * normcdf(x1, 0, sigma) .* (1 - normcdf(x2, 0, sigma)) + W_closed * exp(-0.5 * (x1./sigma).^2) + W_open * exp(-0.5 * (x2./sigma).^2) + P_outlier;

Np = numel(prob_range);
prob = zeros(N, N, Np);
% Apply uncertainty to casted ray
for i = 1:Np
    prob(:,:,i) = double(fun2(prob_range(i) - Zclosed, prob_range(i) - Zopen, sigma));
end
end

