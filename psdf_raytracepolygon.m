function prob = psdf_raytracepolygon( N, M, prob_range, ang, P, Nrays, sigma_x, sigma_y, sigma, P_outlier, W_open, W_closed, W_between )

%% Call self repeatedly if multiple angles are given
if numel(ang) > 1
    prob = psdf_raytracepolygon( N, M, prob_range, ang(1), P, Nrays, sigma_x, sigma_y, sigma, P_outlier, W_open, W_closed, W_between );
    for a = ang(2:end)
        % call self
        p = psdf_raytracepolygon( N, M, prob_range, a, P, Nrays, sigma_x, sigma_y, sigma, P_outlier, W_open, W_closed, W_between );
        prob = prob .* p;
        % Normalize
        prob = bsxfun(@rdivide, prob, max(prob, [], 3));
    end
    return
end

%% Raytrace against geometry
points = raytrace(ang, P, Nrays, sigma_x, sigma_y);

% Compute probability distribution for each ray
prob = psdf_rayprob( N, M, prob_range, ang, points, sigma, P_outlier, W_open, W_closed, W_between );

end
