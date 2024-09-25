function [smoothed_FRs] = guassian_smooth_rasters(rasters, g_sigma, z_score)

smoothed_FRs = NaN(size(rasters));

[n_trials, n_times] = size(rasters);

for t = 1:n_trials
    
    smoothed_FRs(t,:) = gaussian_smooth_data(rasters(t,:), g_sigma);
    
end % of looping over trials

if z_score
    
    smoothed_FRs = zscore(smoothed_FRs, [], 'all');
    
end

end % of function