function [trialinfo, Trodalness, Duration] = get_trialindices(pl2_fname, start_event, end_event, cue_event, choice_event, keep_event)


% open the PL2 file
[OpenedFileName, Version, Freq, Comment,...
    Trodalness, NPW, PreThresh, SpikePeakV,...
    SpikeADResBits, SlowPeakV, SlowADResBits,...
    Duration, DateTime] = plx_information(pl2_fname);

% extract the strobe codes
[nev, evts, evsv] = plx_event_ts(pl2_fname, 257);
% convert strobed words into numbers
% the first event is ALWAYS 9
conversion_amount = evsv(1) - 9;
event_codes = evsv - conversion_amount;
% convert timesteps to milliseconds
event_times = round(evts*1000);


trial_start_indices= find(event_codes == start_event);
trial_end_indices = find(event_codes == end_event);

% the trial start/end indices are repeated in quadlets, let's keep only one
trial_start_indices = trial_start_indices(4:4:end);
trial_end_indices = trial_end_indices(4:4:end);

cue_indices = find(event_codes == cue_event);

choice_indices = find(event_codes == choice_event);
choice_indices = choice_indices(2:2:end);

trial_start_times = event_times(trial_start_indices);
trial_end_times = event_times(trial_end_indices);
cue_times = event_times(cue_indices);
choice_times = event_times(choice_indices);

for t = 1:numel(trial_end_indices)
    
    % get the indices of this trial's starts and ends
    t_start = trial_start_indices(t);
    t_end = trial_end_indices(t);
    
    % get the timestamps of when this trial began and ended
    trialinfo(t,1) = trial_start_times(t);
    trialinfo(t,2) = trial_end_times(t);
    
    % find out whether a choice happened on this trial
    if any(ismember(event_codes(t_start:t_end), keep_event))
        
        % find how many indices away this event was in order to get the
        % time
        
        t_choice_event_ix = choice_indices((choice_indices > t_start) & (choice_indices < t_end));
        t_cue_event_ix = cue_indices((cue_indices > t_start) & (cue_indices < t_end));
        
        % keep track that we want to keep this trial
        trialinfo(t,3) = 1;
        
        % record the timestamp of the choice to ultimately align to
        trialinfo(t,4) = event_times(t_choice_event_ix);
        trialinfo(t,5) = event_times(t_cue_event_ix);
        trialinfo(t,6) = 1;
        
    else
        
        % keep track that we want to keep this trial
        trialinfo(t,3) = 0;
        
        % record the timestamp of the choice to ultimately align to
        trialinfo(t,4) = 0;
        trialinfo(t,5) = 0;
        trialinfo(t,6) = 0;
        
    end % of figuring out which trials to keep
    
end % of cycling over trials for event code aggregation



end % of function