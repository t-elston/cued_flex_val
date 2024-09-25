% high level script to extract spikes from .txt files derived from offline sorter

% needed to make a separate extraction method for King because I sorted his data on the continuous channels in Offline Sorter. 
% I exported the spike times for each unit as a .txt file.

subject = 'King';

spkdir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\Raw data\King\spk_files\';
pl2dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\Raw data\King\tiny_pl2\';
bhvdir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\Raw data\King\bhv2\';
infodir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\Raw data\King\probe_ch_info\';
savedir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\preprocessed data\King\';

spk_files = dir([spkdir '*.txt']);
pl2_files = dir([pl2dir '*.pl2']);
bhv_files = dir([bhvdir '*.bhv2']);
info_files = dir([infodir '*.xlsx']);
[n_spk_files,~] = size(spk_files);

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
    

for f = 1:n_spk_files
    
    fprintf('\n')
    disp(['***loading file: ' spk_files(f).name])
    
    %----------
    % now let's extract the behavioral data
    %----------
    disp(['extracting behavioral data... ', bhv_files(f).name])
    
    % extract behavior and get indices of trial starts, cue times, choice times, and trial ends
    [bhv, ml_align_times] = extractCFVbhv_vRecorded(bhvdir, bhv_files(f).name);
     sname = bhv.fname{1};

    %----------
    % now let's load the sorting notes
    %----------
    disp(['extracting probe details... ', info_files(f).name])
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    probe_notes = readtable([infodir info_files(f).name]); % details of which channels are on which probe
    
    %----------
    % load the tiny PL2 file
    %----------
    disp(['extracting events and LFP... ', pl2_files(f).name])
    pl2_fname = [pl2dir pl2_files(f).name];
   [first_trial_start, Trodalness, Duration] = get_pl2_first_trial_start(pl2_fname, start_event);
    % what was the session length (in milliseconds)?
    session_len = round(Duration*1000);
    
    %----------
    % pull out the spike data
    %----------
    disp(['extracting spike times... ', spk_files(f).name])
    file_spikes = readtable([spkdir spk_files(f).name]);
    file_spikes.Timestamp = round(file_spikes.Timestamp*1000);
    file_spikes.Timestamp(file_spikes.Timestamp==0) = 1;
    ch_nums = unique(file_spikes.Channel_Raw_);
    
    unit_details = unique( file_spikes(:,[1 2]), 'rows');
    n_units = size(unit_details,1);
    disp([num2str(n_units) ' units detected...'])
   
    % convert the monkeylogic timestamps to plexon timestamps
    align_times = ml_align_times + (first_trial_start - ml_align_times(1,1));
    
    % initialize an array that will be the rasters
    n_trials = numel(align_times(:,1));
    %spk_table = zeros(n_trials, offset + offset ,n_units);
    firing_rates=[];
    
    disp('making spike tables...')
    unit_letters = {'a','b','c','d','e','f'};
    unit_names = {};
    unit_channels = [];
    unit_probe_num = [];
    n_units_per_channel = [];
    min_OFC_channel = min(probe_notes.First_ch(contains(probe_notes.Target, 'OFC')));
    max_OFC_channel = max(probe_notes.Last_ch(contains(probe_notes.Target, 'OFC')));
    spk_table=[];
    OFC_ix = [];
    
    spk_table = zeros(n_trials, offset + offset ,n_units);
    
    % now let's loop over the units
    for u = 1:n_units
        
        u_channel = unit_details.Channel_Raw_(u);
        u_num = unit_details.Unit(u);
                
        % find and extract spikes of this unit
        unit_ts = file_spikes.Timestamp((file_spikes.Channel_Raw_ == u_channel) &...
            (file_spikes.Unit == u_num));
        
        % make a cell array of names for each unit
        unit_names{u,1} = [sname, '_', num2str(u_channel), unit_letters{u_num}];
        unit_channels(u) = u_channel;
        
        % which probe was this unit on?
        unit_probe_num(u,1) = find((u_channel >= probe_notes.First_ch) & (u_channel <= probe_notes.Last_ch));
        
        % was this unit in the OFC?
        OFC_ix(u,1) = u_channel <= max_OFC_channel;
        
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
                
                spk_table(t,t_spikes,u) = 1;
                
            end % of filling out spk_table is this trial is used
            
        end % of looping over trials
        
    end % of looping over units
        
    rasters = spk_table;
    rasters(rasters==0) = NaN;
    
    raster_ts = -1*offset: offset-1;
    dense_FR = movmean(spk_table,win_size+1, 2)*1000;
    
    % now let's package all of the relevant info for this session for
    % saving
    trials2keep = bhv.UseTrial == 1;
    
    % now downsample the firing rates
    disp('downsampling spikes...')
    [firing_rates, zFRs, t_mids, mean_FR] = downsample_spikes(dense_FR(trials2keep,:,:), win_size, step_size,...
        raster_ts, n_units);
    
    % compute a drift measure
    trial_means = squeeze(mean(firing_rates, 2));
    prop_valid_trials = sum(trial_means > .5) / sum(trials2keep);
    
    
    % how many units should be saved based on firing rate?
    units2keep = prop_valid_trials > .5; 
    disp(['keeping ' num2str(sum(units2keep)) ' units...']);

    bhv.state(bhv.state==2) = -1;
    bhv.state_type(bhv.state_type==2) = -1;
    
    save_bhv = bhv(trials2keep,3:end);
    bhv_array = table2array(save_bhv);
    bhv_varnames = save_bhv.Properties.VariableNames;
    

if sum(units2keep) > 0
    
    disp('saving data...')
    
    sdata = struct;
    sdata.recfile = sname;
    sdata.probe_notes = probe_notes;
    sdata.unit_names = unit_names(units2keep);
    sdata.unit_probe_nums = unit_probe_num(units2keep);
    sdata.unit_ch_nums = unit_channels(units2keep);
    
    sdata.rasters = rasters(trials2keep,:,units2keep);
    sdata.raw_firing_rates = dense_FR(trials2keep,:,units2keep);
    sdata.raster_ts = raster_ts;
    
    sdata.firing_rates = firing_rates(:,:,units2keep);
    sdata.nrml_FRs = zFRs(:,:,units2keep);
    sdata.t_mids = t_mids;
   
    sdata.OFC_ix = OFC_ix(units2keep);
    
    sdata.bhv = bhv_array;
    sdata.bhv_varnames = bhv_varnames;
    
    % now let's save the data
    save([savedir, sname, '_recdata.mat'],'sdata','-v7.3');
    
else
    
    disp('not enough data to save. ')
    
end
    
    disp('file complete. =]')
    
end % of cycling through files
    