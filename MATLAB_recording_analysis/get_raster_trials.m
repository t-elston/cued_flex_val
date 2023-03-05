function [outrasters] = get_raster_trials(n_trials, rasters, group_ix)

candidate_trials = shuffle(find(group_ix));
try
selected_rasters = rasters(candidate_trials(1:n_trials),:);
catch
    n_trials = sum(group_ix);
    selected_rasters = rasters(candidate_trials(1:n_trials),:);
end
    

selected_rasters(selected_rasters==0) = NaN;

outrasters = [];

for t = 1:n_trials
    
    outrasters(t, :) = selected_rasters(t,:)*t;
end

end