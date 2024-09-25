% high level script to extract spikes from pl2 files and make spike tables
% v04 makes cue-aligned and choice-aligned spike tables

subject = 'Don';

pl2dir = 'E:\CuedFlexVal\sorted_data\';
bhvdir = 'E:\CuedFlexVal\bhv2_files\';
notesdir = 'E:\CuedFlexVal\sorting_notes\';
savedir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\preprocessed data\Don\';

pl2_files= (dir([pl2dir '*.pl2']));
bhv_files = dir([bhvdir '*.bhv2']);
notes_files = dir([notesdir '*.xlsx']);
[n_pl2_files,~] = size(pl2_files);

% some key event codes
keep_event  = 38; % event code indicating a choice was made and the trial should be retained

start_event = 9;
end_event   = 18;
cue_event = 3; % state cues appears on screen
choice_event = 40; % choice options appear on screen
offset = 2000;

% how much downsampling for the spikes do we want?
win_size = 100;
step_size = 25;
    

for f = 1:n_pl2_files
    
    fprintf('\n')
    disp(['***loading file: ' pl2_files(f).name])
    
    %----------
    % now let's extract the behavioral data
    %----------
    disp('extracting behavioral data...')
    
    % extract behavior and get indices of trial starts, cue times, choice times, and trial ends
    [bhv, ml_align_times] = extractCFVbhv_vRecorded(bhvdir, bhv_files(f).name);
    
    
    %----------
    % now let's load the sorting notes
    %----------
    disp('extracting sorting notes...')
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    sorting_notes = readtable([notesdir notes_files(f).name]);
    OFC_ix = double(contains(sorting_notes.Target,'OFC'));
    
    pl2_fname = [pl2dir pl2_files(f).name];
    
    [first_trial_start, Trodalness, Duration] = get_pl2_first_trial_start(pl2_fname, start_event);
    
    % convert the monkeylogic timestamps to plexon timestamps
    align_times = ml_align_times + (first_trial_start - ml_align_times(1,1));
    
    % get some info about the units in this file
    [tscounts, wfcounts, evcounts, slowcounts] = plx_info(pl2_fname,1);
    
    % find the channels with spikes and the number of units on them
    active_channels = find(sum(tscounts) > 0) - 1;
    
    channels_to_check = active_channels;
    
    n_units_per_channel = sum(tscounts(2:end,channels_to_check+1) > 0);
    n_units = sum(n_units_per_channel);
    
    disp([num2str(n_units) ' units detected'])
    
    % what was the session length (in milliseconds)?
    session_len = round(Duration*1000);
    
    % initialize an array that will be the rasters
    n_trials = numel(align_times(:,1));
    spk_table = zeros(n_trials, offset + offset ,n_units);
    firing_rates=[];
    
    disp('making spike tables...')
    
    % now let's loop over channels_to_check
    u_ctr = 0;
    for ch_ix = 1:numel(channels_to_check)
        ch = channels_to_check(ch_ix);
        n_units_on_channel = n_units_per_channel(ch_ix);
        
        for u = 1:n_units_on_channel
            u_ctr = u_ctr+1;
            
            [unit_nts, unit_ts] = plx_ts(pl2_fname, ch, u);
            
            % convert spike timestamps to milliseconds
            unit_ts = round(unit_ts*1000);
            unit_ts(unit_ts==0) = 1;
            
            % now fill out spk_tbl
            for t = 1: n_trials
                
                % check if we use this trial
                if ~isnan(align_times(t,3))
                    
                    t_spikes=[];
                    t_start = align_times(t, 1);
                    t_cue = align_times(t, 2);
                    t_choice = align_times(t, 3);
                    t_end = align_times(t, 4);
                    
                    t_spikes = unit_ts((unit_ts >= t_choice-offset) & ...
                        (unit_ts <= t_choice+offset)) - t_choice + offset;
                    
                    t_spikes(t_spikes==0)=1;
                    
                    spk_table(t,t_spikes,u_ctr) = 1;
                    
                end % of filling out spk_table is this trial is used
                
            end % of looping over trials
            
        end % of looping over units on channel
        
    end % of looping over channels to check
    
    rasters = spk_table;
    rasters(rasters==0) = NaN;
    
    trials2keep = bhv.UseTrial == 1;
    
    raster_ts = -1*offset: offset-1;
    dense_FR = movmean(spk_table,win_size+1, 2)*1000;
    
    % now downsample the firing rates
    disp('downsampling spikes...')
    [firing_rates, zFRs, t_mids, mean_FR] = downsample_spikes(dense_FR(trials2keep,:,:), win_size, step_size,...
        raster_ts, n_units);
    % how many units should be saved based on firing rate?
    units2keep = mean_FR >= 1; % was 2
    disp(['keeping ' num2str(sum(units2keep)) ' units...']);
    
    % now let's package all of the relevant info for this session for
    % saving

    bhv.state(bhv.state==2) = -1;
    bhv.state_type(bhv.state_type==2) = -1;
    
    save_bhv = bhv(trials2keep,3:end);
    bhv_array = table2array(save_bhv);
    bhv_varnames = save_bhv.Properties.VariableNames;
    
    % make a cell array of the names of the units
    sname = bhv.fname{1};
    unit_names = {};
    unit_channels = [];
    unit_probe_num=[];
    for u = 1:numel(units2keep)
        
        unit_names{u,1} = [sname, '_', num2str(sorting_notes.Channel(u)), sorting_notes.Unit{u}];
        unit_probe_num(u,1) = sorting_notes.TowerNumber(u);
        unit_channels(u,1) = sorting_notes.Channel(u);
        
    end % of cycling over units
    
    % make some probe_notes
    probe_nums = unique(sorting_notes.TowerNumber);
    probe_notes = table;
    for p = 1:numel(probe_nums)
        
        if any(contains(sorting_notes.Target(sorting_notes.TowerNumber == probe_nums(p)), 'OFC'))
            
            probe_notes.Target{p} = 'OFC';
        else
            probe_notes.Target{p} = 'HPC';
        end
        
        probe_notes.Tower(p) = probe_nums(p);
        
        probe_notes.First_ch(p) = min(sorting_notes.Channel(sorting_notes.TowerNumber == probe_nums(p)));
        probe_notes.Last_ch(p) = max(sorting_notes.Channel(sorting_notes.TowerNumber == probe_nums(p)));
        
    end % of looping over probes
    
    disp('saving data...')
    
    sdata = struct;
    sdata.recfile = sname;
    sdata.probe_notes = probe_notes;
    sdata.unit_names = unit_names(units2keep);
    sdata.unit_probe_nums = unit_probe_num(units2keep);
    sdata.unit_ch_nums = unit_channels(units2keep);
    
    sdata.rasters = rasters(trials2keep,:,units2keep);
    sdara.raw_firing_rates = dense_FR(trials2keep,:,units2keep);
    sdata.raster_ts = raster_ts;
    
    sdata.firing_rates = firing_rates(:,:,units2keep);
    sdata.nrml_FRs = zFRs(:,:,units2keep);
    sdata.t_mids = t_mids;
    
    sdata.OFC_ix = OFC_ix(units2keep);
    
    sdata.bhv = bhv_array;
    sdata.bhv_varnames = bhv_varnames;
    

    % now let's save the data
    save([savedir, sname, '_recdata.mat'],'sdata','-v7.3');
    
    disp('file complete.')
    
end % of cycling through files
    