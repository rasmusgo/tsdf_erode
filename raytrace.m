function [ points ] = raytrace( ang, P, Nrays, sigma_x, sigma_y )
%% Rotate the polygon
R = [cos(ang) -sin(ang);sin(ang) cos(ang)];
P2 = R * P;

%% Find line segments visible from current view
% Divide the polygon into lines
L_from = P2(:, 1:end-1);
L_to = P2(:, 2:end);

% Flip the lines so they are all pointing to the right, increasing x
L_lo = L_from;
L_hi = L_to;

lr = L_from(1, :) <= L_to(1, :);
L_lo(:, ~lr) = L_to(:, ~lr);
L_hi(:, ~lr) = L_from(:, ~lr);

% Remove minimal segments
valid_idx = L_hi(1,:) - L_lo(1,:) > 1e-10;
L_lo = L_lo(:, valid_idx);
L_hi = L_hi(:, valid_idx);

% Compute the equations for the lines
dp = L_hi - L_lo;
K = dp(2,:) ./ dp(1,:);
M = L_hi(2,:) - K .* L_hi(1,:);

%
[L_leftmost, L_lo_idx] = sort(L_lo(1, :));
[L_rightmost, L_hi_idx] = sort(L_hi(1, :));

% Go through all 'events'
i = 1;
j = 1;
active = false(1, numel(L_lo_idx));
segments = unique([L_leftmost L_rightmost], 'sorted');
L_visible_idx = zeros(size(segments));
for s = 1:numel(segments)
    x = segments(s);
    % Find lines within this portion of the view
    while i <= numel(L_leftmost) && x >= L_leftmost(i)
        active(L_lo_idx(i)) = true;
        i = i + 1;
    end
    while j <= numel(L_rightmost) && x >= L_rightmost(j)
        active(L_hi_idx(j)) = false;
        j = j + 1;
    end
    
    % Determine visible line segment
    y = K(active)*(x + 1e-10) + M(active);

    %find(active)
    if numel(y) > 0
        [y_s, y_idx] = sort(y);
        idx = find(active);
        L_visible_idx(s) = idx(y_idx(1));
    end
end

% Merge segments originating from the same line in the polygon
L_visible_idx = L_visible_idx(1:end-1);
ia = find([true L_visible_idx(1:end-1) ~= L_visible_idx(2:end)]);
u_visible = L_visible_idx(ia);

% NOTE: This code cannot handle reappearing line segments
%[u_visible, ia, ic] = unique(L_visible_idx, 'stable');
%u_visible

u_segments = [segments(ia) segments(end)];

L_visible_lo = [u_segments(1:end-1); K(u_visible) .* u_segments(1:end-1) + M(u_visible)];
L_visible_hi = [u_segments(2:end);   K(u_visible) .* u_segments(2:end)   + M(u_visible)];

% Remove minimal segments
valid_idx = L_visible_hi(1,:) - L_visible_lo(1,:) > 1e-10;
L_visible_lo = L_visible_lo(:,valid_idx);
L_visible_hi = L_visible_hi(:,valid_idx);
K = K(u_visible(valid_idx));
M = M(u_visible(valid_idx));

% Sample x-values to trace
X_ideal = linspace(-10, 10, Nrays);
X = X_ideal + sigma_x * randn(1, Nrays);

% Find intersection points of those rays
Y = nan(1, Nrays);
j = 1;
jmax = size(L_visible_hi, 2);
for i=1:numel(X)
    while j <= jmax && X(i) > L_visible_hi(1,j)
        j = j + 1;
    end
    if j > jmax
        break
    end
    if X(i) > L_visible_lo(1,j)
        Y(i) = K(j) * X(i) + M(j);
    end
end

% Add noise to Y
Y = Y + sigma_y * randn(1, Nrays);

% Undo rotation and collect final result
points = R' * [X_ideal ; Y];

end
