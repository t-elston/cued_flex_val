% high level script to extract spikes from pl2 files and make spike tables
% v02 makes cue-aligned and choice-aligned spike tables

pl2dir = 'E:\CuedFlexVal\sorted_data\'; 
bhvdir = 'E:\CuedFlexVal\bhv2_files\'; 
notesdir = 'E:\CuedFlexVal\sorting_notes\';
savedir = 'C:\Users\Thomas Elston\Documents\PYTHON\CuedFlexVal\recdata\';

pl2_files= dir([pl2dir '*.pl2']);
bhv_files = dir([bhvdir '*.bhv2']);
notes_files = dir([notesdir '*.xlsx']);
[n_pl2_files,~] = size(pl2_files);


% some key event codes
keep_event  = 38; % event code indicating a choice was made and the trial should be retained

start_event = 9;
end_event   = 18;
cue_event = 3; % state cues appears on screen
choice_event = 40; % choice options appear on screen
cue_psth_pos_offset = 500;
cue_psth_neg_offset = 500;
choice_psth_pos_offset = 1000;
choice_psth_neg_offset = 1000;


% how much downsampling for the spikes do we want?
win_size = 100;
step_size = 10;

% make filter for bandpassing LFP into theta range
theta_filt = designfilt('bandpassfir', ... response type
    'FilterOrder',1000, ... filter order
    'CutoffFrequency1', 4,... lower range
    'CutoffFrequency2', 8,... upper range
    'SampleRate', 1000); ... sampling rate


for f = 1:n_pl2_files
    
    pl2_fname = [pl2dir pl2_files(f).name];
    
    fprintf('\n')  
    disp(['***loading file: ' pl2_files(f).name]) 
    
    % get indices of trial starts, ends, cue times, and choice times in
    % milliseconds (absolute recording time)
    [trialinfo, Trodalness, Duration] = get_trialindices(pl2_fname, start_event, end_event, cue_event, choice_event, keep_event);

    % now I want to make spike tables
    % get some info about the units in this file
    [tscounts, wfcounts, evcounts, slowcounts] = plx_info(pl2_fname,1);
    
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
    cue_spk_table = NaN(n_trials,cue_psth_neg_offset + cue_psth_pos_offset,total_n_units);
    choice_spk_table = NaN(n_trials,choice_psth_neg_offset + choice_psth_pos_offset,total_n_units);

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
            
            [unit_nts, unit_ts] = plx_ts(pl2_fname, ch, u);
            
            % convert spike timestamps to milliseconds
            unit_ts = round(unit_ts*1000);
            unit_ts(unit_ts==0) = 1;
            
            % get the ISIs for each unit
            u_ISI{u_ctr,1} = diff(unit_ts);
            u_mean_ISIs(u_ctr,1) = mean(diff(unit_ts));
            
            % fill out a blank array with 1s and zeros where spikes occured
            spk_array = zeros(1,session_len);
            spk_array(unit_ts) = 1;
            
            % now derive firing rates
            spk_firing_rates = movmean(spk_array,win_size+1)*1000;
            
            % now fill out spk_tbl
            for t = 1: n_trials
                    
                % check if we use this trial
                if trialinfo(t,3) == 1
                    
                    t_cue_time = trialinfo(t, 5);
                    t_choice_time = trialinfo(t, 4);

                    choice_spk_table(t,:,u_ctr) = spk_array(t_choice_time - choice_psth_neg_offset : t_choice_time + choice_psth_pos_offset-1);
                    cue_spk_table(t,:,u_ctr) = spk_array(t_cue_time - cue_psth_neg_offset : t_cue_time + cue_psth_pos_offset-1);
                    
                    choice_raw_firing_rates(t,:,u_ctr) = spk_firing_rates(t_choice_time - choice_psth_neg_offset : t_choice_time + choice_psth_pos_offset-1);
                    cue_raw_firing_rates(t,:,u_ctr) = spk_firing_rates(t_cue_time - cue_psth_neg_offset : t_cue_time + cue_psth_pos_offset-1);
                    
                end % of filling out spk_table is this trial is used

            end % of looping over trials

        end % of looping over units on channel
       
    end % of looping over channels to check
    
    cue_raw_timesteps = -1*cue_psth_neg_offset: cue_psth_pos_offset-1;
    choice_raw_timesteps = -1*choice_psth_neg_offset: choice_psth_pos_offset-1;

    
    % now downsample the firing rates
    disp('downsampling spikes...')
    [cue_down_FRs, cue_zFRs, cue_down_ts, cue_mean_FR] = downsample_spikes(cue_raw_firing_rates, win_size, step_size, cue_raw_timesteps, total_n_units); 
    [choice_down_FRs, choice_zFRs, choice_down_ts, choice_mean_FR] = downsample_spikes(choice_raw_firing_rates, win_size, step_size, choice_raw_timesteps, total_n_units); 
    
    % how many units should be saved based on firing rate?
    units2keep = choice_mean_FR > 1 | cue_mean_FR > 1;
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
    OFC_ix = double(contains(sorting_notes.Target,'OFC'));

    % now let's package all of the relevant info for this session for
    % saving
    
    trials2keep = bhv.UseTrial == 1;
    
    bhv.state(bhv.state==2) = -1;
    bhv.state_type(bhv.state_type==2) = -1;

    save_bhv = bhv(trials2keep,3:end);
    bhv_array = table2array(save_bhv);
    bhv_varnames = save_bhv.Properties.VariableNames;
    
    % make a cell array of the names of the units
      sname = bhv.fname{1};
      unit_names = {};
      unit_channels = [];
      for u = 1:numel(units2keep)
          
          unit_names{u,1} = [sname, '_', num2str(sorting_notes.Channel(u)), sorting_notes.Unit{u}];
          unit_channels(u,1) = sorting_notes.Channel(u);
          
      end % of cycling over units
      
      
      disp('extracting and bandpassing LFP...')
      % bandpass LFP into theta and then extract phase and amplitude
      lfp_channels_to_check = channels_to_check(find(n_units_per_channel > 0));
      % extract the LFP and find channels with good signal

    % get details about the LFP channel names
    [~,file_ch_names] = plx_adchan_names(pl2_fname);
    file_ch_names = cellstr(file_ch_names);
    file_ch_names = file_ch_names(cellfun(@(x) ~isempty(strfind(x,'FP')),file_ch_names)); % only keep FP channels
    
    cue_theta_amp = NaN(n_trials,cue_psth_neg_offset + cue_psth_pos_offset,numel(lfp_channels_to_check));
    cue_theta_phases = NaN(n_trials,cue_psth_neg_offset + cue_psth_pos_offset,numel(lfp_channels_to_check));
    choice_theta_amp = NaN(n_trials,choice_psth_neg_offset + choice_psth_pos_offset,numel(lfp_channels_to_check));
    choice_theta_phases = NaN(n_trials,choice_psth_neg_offset + choice_psth_pos_offset,numel(lfp_channels_to_check));
    
    cue_theta_wave = NaN(n_trials,cue_psth_neg_offset + cue_psth_pos_offset,numel(lfp_channels_to_check));
    choice_theta_wave = NaN(n_trials,choice_psth_neg_offset + choice_psth_pos_offset,numel(lfp_channels_to_check));

    
    lfp_in_ofc = [];
    lfp_ch_num = [];

    for ch_ix = 1:numel(lfp_channels_to_check)
        
        lfp_in_ofc(ch_ix,1) = double(contains(sorting_notes.Target(min(find(sorting_notes.Channel == lfp_channels_to_check(ch_ix)))),'OFC'));
        lfp_ch_num(ch_ix,1) = lfp_channels_to_check(ch_ix);
        
        ChLFP = PL2Ad(pl2_fname, file_ch_names{lfp_channels_to_check(ch_ix)});
        
        theta_signal = fftfilt(theta_filt.Coefficients,ChLFP.Values);
        
        theta_amp = abs(hilbert(theta_signal));
        theta_phase = angle(hilbert(theta_signal));
        
        % adjust time to account for filter lag
        LFPtime = 1:numel(ChLFP.Values); %ms
        lag = 1+(1+1)*ChLFP.ADFreq/2;
        
        adjusted_theta_wave = theta_signal(lag:end);
        adjusted_theta_amp = theta_amp(lag:end);
        adjusted_theta_phase = theta_phase(lag:end);
        
        % now chop the adjusted phase and amplitudes into trials
        for t = 1: n_trials
            % check if we use this trial
            
            if trialinfo(t,3) == 1
                
                t_cue_time = trialinfo(t, 5);
                t_choice_time = trialinfo(t, 4);
                
                cue_theta_amp(t,:,ch_ix) = adjusted_theta_amp(t_cue_time - cue_psth_neg_offset : t_cue_time + cue_psth_pos_offset-1);
                cue_theta_phase(t,:,ch_ix) = adjusted_theta_phase(t_cue_time - cue_psth_neg_offset : t_cue_time + cue_psth_pos_offset-1);
                
                choice_theta_amp(t,:,ch_ix) = adjusted_theta_amp(t_choice_time - choice_psth_neg_offset : t_choice_time + choice_psth_pos_offset-1);
                choice_theta_phase(t,:,ch_ix) = adjusted_theta_phase(t_choice_time - choice_psth_neg_offset : t_choice_time + choice_psth_pos_offset-1);
                
                cue_theta_wave(t,:,ch_ix) = adjusted_theta_wave(t_cue_time - cue_psth_neg_offset : t_cue_time + cue_psth_pos_offset-1);
                choice_theta_wave(t,:,ch_ix) = adjusted_theta_wave(t_choice_time - choice_psth_neg_offset : t_choice_time + choice_psth_pos_offset-1);
                
            end
        end
    end
    
    % now look at phase-locking
    disp('assessing unit phase locking...')
    [~,cue_win_start] = min(abs(cue_raw_timesteps - -100));
    [~,cue_win_end] = min(abs(cue_raw_timesteps - 600));
    [~,choice_win_start] = min(abs(choice_raw_timesteps - -100));
    [~,choice_win_end] = min(abs(choice_raw_timesteps - 600));
    
    [cue_phase_locking] = CFV_unit_phase_locking(cue_spk_table(trials2keep,cue_win_start:cue_win_end,:),...
                                            cue_theta_phase(trials2keep,cue_win_start:cue_win_end,:),... 
                                            cue_theta_wave(trials2keep,cue_win_start:cue_win_end,:),unit_channels, lfp_ch_num);
      
    [choice_phase_locking] = CFV_unit_phase_locking(choice_spk_table(trials2keep,choice_win_start:choice_win_end,:),...
                                            choice_theta_phase(trials2keep,choice_win_start:choice_win_end,:),... 
                                            choice_theta_wave(trials2keep,choice_win_start:choice_win_end,:),unit_channels, lfp_ch_num);                                    
    
                                        
                                        
    % clear a bunch of variables  
    cue_theta_wave=[];
    cue_theta_amp=[];
    cue_theta_phase=[];
    choice_theta_wave=[];
    choice_theta_amp=[];
    choice_theta_phase=[];
    ChLFP = [];
                                        
    disp('saving data...')
    
    
    sdata = struct; 
    sdata.recfile = sname;
    sdata.unit_names = unit_names(units2keep);
    sdata.unit_ch_nums = unit_channels(units2keep);
    sdata.choice_rasters = choice_spk_table(trials2keep,:,units2keep);
    sdata.cue_rasters = cue_spk_table(trials2keep,:,units2keep);
      
    sdata.raw_choice_FRs = choice_raw_firing_rates(trials2keep,:,units2keep);
    sdata.raw_cue_FRs = cue_raw_firing_rates(trials2keep,:,units2keep);

    sdata.down_choice_FRs = choice_down_FRs(trials2keep,:,units2keep);
    sdata.down_cue_FRs = cue_down_FRs(trials2keep,:,units2keep);

    sdata.cue_raw_times = cue_raw_timesteps;
    sdata.choice_raw_times = choice_raw_timesteps;

    sdata.cue_down_times = cue_down_ts;
    sdata.choice_down_times = choice_down_ts;

    sdata.choice_zFRs = choice_zFRs(trials2keep,:,units2keep);
    sdata.cue_zFRs = cue_zFRs(trials2keep,:,units2keep);
 
    sdata.unit_ISIs = u_ISI(units2keep);
    sdata.mean_ISIs = u_mean_ISIs(units2keep);
    
    sdata.phase_locking_varnames = {'pval','mean_angle','vector_len','zval'};
    sdata.cue_phase_locking = cue_phase_locking(units2keep,:);
    sdata.choice_phase_locking = choice_phase_locking(units2keep,:);

    sdata.OFC_ix = OFC_ix(units2keep);

    sdata.bhv = bhv_array;
    sdata.bhv_varnames = bhv_varnames;
    
%     sdata.cue_theta_amp = cue_theta_amp;
%     sdata.cue_theta_phase = cue_theta_phase;
%     sdata.choice_theta_amp = choice_theta_amp;
%     sdata.choice_theta_phase = choice_theta_phase;
%     sdata.lfp_in_ofc = lfp_in_ofc;
%     sdata.lfp_ch_num = lfp_ch_num;
%     
    
    
    % now let's save the data
    save([savedir, sname, '_recdata.mat'],'-struct','sdata','-v7.3');
    
    disp('file complete.')

end % of cycling through files
    