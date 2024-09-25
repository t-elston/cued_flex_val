function [pl2_event_times] = get_pl2_trial_event_timestamps(pl2_file, event_ids, n_kept_trials)

pl2_event_times = table;

% extract the strobe codes
[nev, evts, evsv] = plx_event_ts(pl2_file, 257);

% convert strobed words into numbers
% the first event is ALWAYS 9
conversion_amount = evsv(1) - 9;
event_codes = evsv - conversion_amount;
% convert timesteps to milliseconds
event_times = round(evts*1000);

t_starts = strfind(event_codes',[9, 9, 9, 9]);
t_ends = strfind(event_codes',[18, 18, 18, 18]);

% now loop over all trials
for t = 1:numel(t_starts)
    
    trial_events = event_codes(t_starts(t) : t_ends(t));
    trial_event_times = event_times(t_starts(t) : t_ends(t));
    
    choice_made(t,1) = ~isempty(find(trial_events == event_ids.keep));
    
    if choice_made(t)
    
        cue_on_ix  = find(trial_events == event_ids.cue_on);
        pics_on_ix = min(find(trial_events == event_ids.pics_on));
        
        cue_on(t,1) = trial_event_times(cue_on_ix);
        pics_on(t,1) = trial_event_times(pics_on_ix);
        
    end

end % of looping over trials 

% only save the trials where a choice was actually made
pl2_event_times.cue_on = cue_on(choice_made);
pl2_event_times.pics_on = pics_on(choice_made);


end % of function