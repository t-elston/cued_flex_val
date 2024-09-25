function [lfp_trials] = chop_lfp_into_trials(lfp_data, event_times, offset)

n_trials = numel(event_times.pics_on);

lfp_trials = NaN(n_trials, 2*offset);

% loop over individual trials
for t = 1:numel(event_times.pics_on)
    
    lfp_trials(t, :) = lfp_data(event_times.pics_on(t) - offset : event_times.pics_on(t) + offset-1);

end % of looping over trials

end % of function