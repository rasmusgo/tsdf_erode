function logodds = occu_raytracepolygon( N, M, ang, P, Nrays, sigma_x, sigma_y, sigma, P_outlier )

%% Call self repeatedly if multiple angles are given
if numel(ang) > 1
    prior = log(ones(N));
    logodds = occu_raytracepolygon( N, M, ang(1), P, Nrays, sigma_x, sigma_y, sigma, P_outlier );
    for a = ang(2:end)
        % call self
        p = occu_raytracepolygon( N, M, a, P, Nrays, sigma_x, sigma_y, sigma, P_outlier );
        logodds = logodds + p - prior;
    end
    return
end

%% Raytrace against geometry
points = raytrace(ang, P, Nrays, sigma_x, sigma_y);

% Compute probability distribution for each ray
logodds = occu_rayprob( N, M, ang, points, sigma, P_outlier );

end
