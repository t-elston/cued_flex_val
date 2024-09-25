% compute_phase_locking.m

% This script computes the within and across-region phase-locking units in a recording.
% The relevant data are then added to the struct (sdata) containing all data for a given session and re-saved. 

rec_dir = 'C:/Users/Thomas Elston/Documents/MATLAB/Projects/CuedFlexVal/preprocessed data/';
raw_dir = 'D:/Raw data/';

% get subject-level directories 
subject_folders = dir(rec_dir);
subject_folders = subject_folders(~ismember({subject_folders(:).name},{'.','..'}));
subject_folders = subject_folders(cell2mat({subject_folders(:).isdir}));

fprintf('\nrunning LFP phase and coherence analysis');

% now loop over each subject 
for s = 1:numel(subject_folders)
    
    % now find the preprocessed files for this subject
    rec_fnames = dir([rec_dir '/' subject_folders(s).name '/units/']);
    rec_fnames = rec_fnames(~ismember({rec_fnames(:).name},{'.','..'}));
    
    % find the lfp_dir associated with this subject
    lfp_fnames = dir([raw_dir '/' subject_folders(s).name]);
    lfp_fnames = lfp_fnames(~ismember({lfp_fnames(:).name},{'.','..'}));
    
    % now loop over each recording
    for r = 1:numel(rec_fnames)
        
        fprintf(['\nprocessing ' rec_fnames(r).name '\n']);

        % load this recording file
        load([rec_fnames(r).folder '/' rec_fnames(r).name]);

        % extract behavioral data
        bhv = array2table(sdata.bhv, 'VariableNames', sdata.bhv_varnames);
        bhv_varnames = sdata.bhv_varnames;
        
        n_kept_trials = numel(sdata.bhv(:,1));
        
        % find the pl2 associated with this file
        target_pl2 = find(contains({lfp_fnames.name}, sdata.recfile));
        
        % extract the relevant trials and timestamps for the cue on and choice on
        event_ids = table; 
        event_ids.t_start = 9; 
        event_ids.cue_on = 3; 
        event_ids.pics_on = 40;
        event_ids.keep = 38;
        event_ids.t_end = 18; 
        
        pl2_fname = [lfp_fnames(target_pl2).folder '/' lfp_fnames(target_pl2).name];
        [pl2_event_times] = get_pl2_trial_event_timestamps(pl2_fname, event_ids, n_kept_trials);
        
        % get details about the LFP channel names
        [~,pl2_ch_names] = plx_adchan_names(pl2_fname);
        pl2_ch_names = cellstr(pl2_ch_names);
        
        pl2_ch_names = pl2_ch_names(cellfun(@(x) ~isempty(strfind(x,'FP')),pl2_ch_names)); % only keep FP channels
         
        fprintf('assessing LFP phase alignment...\n');
        FOI = [ 3, 8;
               10, 15;
               15, 30;
               40, 70];
        
        % get trialwise lfp phase for the different FQ bands. Also get each channels PLV
        [LFP_phase, ch_details] = get_LFP_phase(sdata, pl2_fname, pl2_ch_names, pl2_event_times, FOI);
        
        % now check and see which channels have theta
        fprintf('computing PSDs and identifying channels with theta...\n');
        ch_details = compute_channel_PSDs_v3(ch_details, pl2_fname);
        
        % for the coherence analysis, only use one channel per probe - the one with the biggest theta signal
        fprintf('identifying highest S/N LFP channels...\n');
        probe_ids = unique(ch_details.Tower);
        channels2use = table;
        for p = 1:numel(probe_ids)
            
            this_probe_data = ch_details(ch_details.Tower == probe_ids(p),:);
            [~,max_theta_ix] = max(this_probe_data.theta_power);
            channels2use(p,:) = this_probe_data(max_theta_ix,:);
            
        end % of looping over probes
        
    % loop over channel pairs in this session and compute relevant coherograms
    if ~isempty(channels2use)
                
        [trial_coh, z_coh, coh_ts, PDC, Granger] = compute_coherograms(channels2use, pl2_fname, pl2_event_times);

    end
     
        % save the data
        fprintf('saving data...\n');
        s_name = [rec_dir subject_folders(s).name '/LFP/' sdata.recfile '_LFP'];
        save(s_name, 'LFP_phase','ch_details','trial_coh','z_coh','coh_ts','PDC',...
            'Granger','bhv','bhv_varnames', '-v7.3');
        
        LFP_phase = [];
        
    end % of looping over recordings
    
end % of looping over subjects

fprintf('\n all done :] \n ');