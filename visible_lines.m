function [ L ] = visible_lines( ang, P )
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
valid_idx = L_hi(1,:) - L_lo(1,:) > 1e-12;
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
    
    L_active_lo = L_lo(:, active);
    L_active_hi = L_hi(:, active);
    
    % Determine visible line segment
    dp = L_active_hi - L_active_lo;
    k = dp(2,:) ./ dp(1,:);
    m = L_active_hi(2,:) - k .* L_active_hi(1,:);
    y = k*(x+0.001) + m;

    y = K(active)*(x+0.00001) + M(active);

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

% Undo rotation and collect final result
L = [R' * L_visible_lo ; R' * L_visible_hi];

% plot([L_visible_lo(1,:); L_visible_hi(1,:)], [L_visible_lo(2,:); L_visible_hi(2,:)])
% axis equal
% axis([-10,10, -10,10])


end
