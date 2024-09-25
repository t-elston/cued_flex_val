function [phase_alignment] = compute_LFP_phase_alignment_v2(inwaves)
% inputs:
% - inwaves = n_trials x n_times array comprised of LFP phases

% outputs
% - phase_alignment = 1D n_times long array with mean phase alignment for each timestep

[n_trials, n_times] = size(inwaves);

% compute phase alignment at each timepoint
for t = 1:n_times
    phase_alignment(1, t) = abs(sum(exp(i*((inwaves(:,t))))))/n_trials;
end % of looping over timesteps

end % of function