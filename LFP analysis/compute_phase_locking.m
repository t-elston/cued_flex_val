% compute_phase_locking.m

% This script computes the within and across-region phase-locking units in a recording.
% The relevant data are then added to the struct (sdata) containing all data for a given session and re-saved. 

rec_dir = 'C:/Users/Thomas Elston/Documents/MATLAB/Projects/CuedFlexVal/preprocessed data/';
raw_dir = 'C:/Users/Thomas Elston/Documents/MATLAB/Projects/CuedFlexVal/raw data/';

% get subject-level directories 
subject_folders = dir(rec_dir);
subject_folders = subject_folders(~ismember({subject_folders(:).name},{'.','..'}));
subject_folders = subject_folders(cell2mat({subject_folders(:).isdir}));

fprintf('\nrunning phase-locking analysis');

% now loop over each subject 
for s = 1:numel(subject_folders)
    
    % now find the preprocessed files for this subject
    rec_fnames = dir([rec_dir '/' subject_folders(s).name '/units/']);
    rec_fnames = rec_fnames(~ismember({rec_fnames(:).name},{'.','..'}));
    
    % find the lfp_dir associated with this subject
    lfp_fnames = dir([raw_dir '/' subject_folders(s).name '/tiny_pl2']);
    lfp_fnames = lfp_fnames(~ismember({lfp_fnames(:).name},{'.','..'}));
    
    % now loop over each recording
    for r = 1:numel(rec_fnames)
        
        fprintf(['\nprocessing ' rec_fnames(r).name '\n']);

        % load this recording file
        load([rec_fnames(r).folder '/' rec_fnames(r).name]);
        
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
          
        
        fprintf('doing phase-locking across multiple fq bands...\n');
        % now actually do the phase-locking analysis
        FOI = [4, 8;
               8, 12;
               15, 30;
               40, 70];
               
        [cue_pvals, cue_PLV, cue_mean_angle, ...
          pics_pvals, pics_PLV, pics_mean_angle] = do_phase_locking_analysis(sdata, pl2_fname, pl2_ch_names, pl2_event_times, FOI);
        
        
         %--------------------------------------------
         %   add information to sdata to be resaved
         %--------------------------------------------
         fprintf('re-saving the data.\n');
         
         sdata.pl2_event_times = pl2_event_times;
         
         sdata.phase_locking_bands = {'theta','alpha','beta','gamma'};
         sdata.phase_locking_freqs = FOI;
         sdata.cue_phase_locking = cue_pvals; 
         sdata.PL_cue_mean_angle = cue_mean_angle; 
         sdata.PLV_cue = cue_PLV; 
         
         sdata.choice_phase_locking = pics_pvals; 
         sdata.PL_pics_mean_angle = pics_mean_angle; 
         sdata.PLV_pics = pics_PLV; 
        
        % re-save the data
        save([rec_fnames(r).folder '/' rec_fnames(r).name],'sdata','-v7.3');
        
    end % of looping over recordings
    
end % of looping over subjects

fprintf('\n all done :] \n ');