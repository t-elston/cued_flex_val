% high level script to extract spikes from pl2 files and make spike tables

pl2dir = 'E:\CuedFlexVal\sorted_data\'; 
bhvdir = 'E:\CuedFlexVal\bhv2_files\'; 
notesdir = 'E:\CuedFlexVal\sorting_notes\';
savedir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\recording_data\';

pl2_files= dir([pl2dir '*.pl2']);
bhv_files = dir([bhvdir '*.bhv2']);
notes_files = dir([notesdir '*.xlsx']);
[n_pl2_files,~] = size(pl2_files);


% some key event codes
keep_event  = 38; % event code indicating a choice was made and the trial should be retained
start_event = 9;
end_event   = 18;
align_event = 40; % choice options appear on screen
cue_on = 3;
psth_pos_offset = 2000;
psth_neg_offset = 1999;

% how much downsampling do we want?
win_size = 100;
step = 10;


for f = 1:n_pl2_files
    
    pl2_fname = [pl2dir pl2_files(f).name];
    
    fprintf('\n')  
    disp(['***loading file: ' pl2_files(f).name])
    
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
    
    % find the trial start/end times and whether to include the trial
    trialinfo = [];
    
    trial_start_indices= find(event_codes == start_event);
    trial_end_indices = find(event_codes == end_event);
    
    % the trial start/end indices are repeated in quadlets, let's keep only one
    trial_start_indices = trial_start_indices(4:4:end);
    trial_end_indices = trial_end_indices(4:4:end);
    
    align_indices = find(event_codes == align_event);
    align_indices = align_indices(2:2:end);
    
    
    trial_start_times = event_times(trial_start_indices);
    trial_end_times = event_times(trial_end_indices);
    align_times = event_times(align_indices);
    cue_indices = find(event_codes == cue_on);
        
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
                
                t_align_event_ix = align_indices((align_indices > t_start) & (align_indices < t_end));
                t_cue_event_ix = cue_indices((cue_indices > t_start) & (cue_indices < t_end));
                
                % keep track that we want to keep this trial
                trialinfo(t,3) = 1;
                
                % record the timestamp of the choice to ultimately align to
                trialinfo(t,4) = event_times(t_align_event_ix); 
                trialinfo(t,5) = event_times(t_cue_event_ix); 
                
            else
                
                % keep track that we want to keep this trial
                trialinfo(t,3) = 0;
                
                % record the timestamp of the choice to ultimately align to
                trialinfo(t,4) = 0; 

            end % of figuring out which trials to keep
               
    end % of cycling over trials for event code aggregation
    
    
    % now I want to make spike tables
    % get some info about the units in this file
    [tscounts, wfcounts, evcounts, slowcounts] = plx_info(OpenedFileName,1);
    
    % find the channels with spikes and the number of units on them
    active_channels = find(sum(tscounts) > 0) - 1;
    
    % take trodality into account - if it's stereotrode, only use the odd
    % numbered channels
    if numel(Trodalness) == 2
        channels_to_check = active_channels(logical(mod(active_channels, numel(Trodalness))));
    else
        channels_to_check = active_channels;
    end
        
    n_units_per_channel = sum(tscounts(2:end,channels_to_check+1) > 0);
    total_n_units = sum(n_units_per_channel);
    
    disp([num2str(total_n_units) ' units detected'])
    
    % what was the session length (in milliseconds)?
    session_len = round(Duration*1000);
    
    % initialize an array that will be the rasters
    n_trials = numel(trialinfo(:,1));
    spk_table = NaN(n_trials,psth_neg_offset + psth_pos_offset,total_n_units);
    
    % make a single array to loop over for extracting units off from
    % channels
    
    disp('making spike tables...')
    
    % now let's loop over channels_to_check
    u_ISI={};
    u_mean_ISIs=[];
    u_ctr = 0;
    for ch_ix = 1:numel(channels_to_check)
        ch = channels_to_check(ch_ix);
        n_units_on_channel = n_units_per_channel(ch_ix);
        
        for u = 1:n_units_on_channel
             u_ctr = u_ctr+1;
            
            [unit_nts, unit_ts] = plx_ts(OpenedFileName, ch, u);
            
            % convert spike timestamps to milliseconds
            unit_ts = round(unit_ts*1000);
            unit_ts(unit_ts==0) = 1;
            
            % get the ISIs for each unit
            u_ISI{u_ctr,1} = diff(unit_ts);
            u_mean_ISIs(u_ctr,1) = mean(diff(unit_ts));
            
            
            
            % fill out a blank array with 1s and zeros where spikes occured
            spk_array = zeros(1,session_len);
            spk_array(unit_ts) = 1;
            
            % now fill out spk_tbl
            for t = 1: n_trials
                    
                % check if we use this trial
                if trialinfo(t,3) == 1
                    
                    t_choice_time = trialinfo(t, 4);
                    
                    spk_table(t,:,u_ctr) = spk_array(t_choice_time - psth_neg_offset : t_choice_time + psth_pos_offset-1);
                    
                end % of filling out spk_table is this trial is used

            end % of looping over trials

        end % of looping over units on channel
       
    end % of looping over channels to check
    
    
    % OK, now convert spk_table to firing rates
    
    raw_firing_rates = movmean(spk_table,100,2) * 1000;
    raw_timesteps = -1*psth_neg_offset: psth_pos_offset;
    
    % now downsample the firing rates
    % make array of window start/ends for downsampling
    win_centers = win_size/2:25:numel(spk_table(1,:,1));
    win_starts = win_centers-win_size/2;
    win_starts(win_starts==0)=1;
    win_ends   = win_centers+win_size/2;
    win_ends(win_ends > numel(raw_firing_rates(1,:,1))) = numel(raw_firing_rates(1,:,1));
    
    down_firing_rates = [];
    down_timesteps = [];
    z_firing_rates=[];
    
    mean_firing_rates = [];
    for u = 1:total_n_units
        
        mean_firing_rates(u,1) = nanmean(raw_firing_rates(:,:,u),'all');
        
        for w = 1:numel(win_centers)
            down_firing_rates(:,w,u) =  nanmean(raw_firing_rates(:,win_starts(w) : win_ends(w),u),2);
            
            down_timesteps(w) = nanmean(raw_timesteps(win_starts(w) : win_ends(w)));
        end % of cylcing over win_centers
        
        % zscore the firing rates
        u_mean = nanmean(down_firing_rates(:,:,u),'all');
        u_std  = nanstd(down_firing_rates(:,:,u),[],'all');
        
        z_firing_rates(:,:,u) = (down_firing_rates(:,:,u) - u_mean) / u_std;
        
    end % of downsampling firing rates
    

    % how many units should be saved based on firing rate?
      units2keep = mean_firing_rates >= 1 & u_mean_ISIs <= 500;
      disp(['keeping ' num2str(sum(units2keep)) ' units...']);
      

    %----------
    % now let's extract the behavioral data
    %----------
    disp('extracting behavioral data...')

    [bhv, ml_align_times] = extractCFVbhv_vRecorded(bhvdir, bhv_files(f).name);
    
    %----------
    % now let's load the sorting notes
    %----------
    disp('extracting sorting notes...')
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    sorting_notes = readtable([notesdir notes_files(f).name]);
     

    
    % now let's package all of the relevant info for this session for
    % saving
    
    disp('saving data...')
    
    sname = bhv.fname{1};
    
    sdata = struct; 
    sdata.recfile = sname;
    sdata.rasters = spk_table(:,:,units2keep);
    sdata.raw_firing_rate = raw_firing_rates(:,:,units2keep);
    sdata.raw_times = raw_timesteps;
    sdata.unit_ISIs = u_ISI;
    sdata.mean_ISIs = u_mean_ISIs;
    sdata.down_firing_rates = down_firing_rates(:,:,units2keep);
    sdata.down_timesteps = down_timesteps;
    sdata.z_firing_rates = z_firing_rates;
    sdata.bhv = bhv;
    sdata.sorting_notes = sorting_notes(units2keep,:);
    
    % now let's save the data
    save([savedir, sname, '_recdata.mat'],'-struct','sdata','-v7.3');
    
    disp('file complete.')

end % of cycling through files
    