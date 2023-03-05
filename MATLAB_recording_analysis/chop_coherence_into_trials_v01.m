function [trial_data, iti_data, coh_ts] = chop_coherence_into_trials_v01(coh,ts,trialinfo)


% what's the new sampling frequency of the coherence?
coh_fs = round(mean(diff(ts)));

neg_offset = 2000;
pos_offset = 2000;

% pull out coherence data aligned to pics on
trialinfo(trialinfo(:,3)==0,:) = [];   
choice_times = trialinfo(:,4);
trial_starts = trialinfo(:,1);
n_trials = numel(choice_times);

coh_ts = neg_offset*-1 : coh_fs : pos_offset;


% cycle through the trials
for t = 1:n_trials
          
        % find the closest trial start/stop times in the spectrogram timestamps
        [~,ts_start] = min(abs(ts - (choice_times(t)-neg_offset)));
        [~,ts_end] = min(abs(ts - (choice_times(t)+pos_offset))); 
                
        % find the closest ITI start/stop times in the spectrogram timestamps
        [~,iti_start] = min(abs(ts - (trial_starts(t)-500)));
        [~,iti_end] = min(abs(ts - (trial_starts(t)-200)));
        
        trial_data(:,:,t) = coh(:,ts_start: ts_end);
        iti_data(:,:,t) = coh(:,iti_start: iti_end);% had been +15
        
end % of cycling through trials


end % of function