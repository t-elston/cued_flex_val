% unit_selectivity_analysis.m

%-----------------------------
% This function looks to see whether single neurons encoded context and other variables during the cue
% and choice phases of the task.
%
% The function also saves trial-averaged data for each neuron
%-----------------------------
rec_dir = 'C:/Users/Thomas Elston/Documents/MATLAB/Projects/CuedFlexVal/preprocessed data/';
raw_dir = 'C:/Users/Thomas Elston/Documents/MATLAB/Projects/CuedFlexVal/raw data/';

% get subject-level directories
subject_folders = dir(rec_dir);
subject_folders = subject_folders(~ismember({subject_folders(:).name},{'.','..'}));
subject_folders = subject_folders(cell2mat({subject_folders(:).isdir}));

fprintf('\nrunning unit selectivity analysis \n');

% initialize a table to accumulate all data into
selectivity_table = table;
unit_data = struct;
u_ctr = 0;

% loop over each subject
for s = 1:numel(subject_folders)
    
    % find the preprocessed files for this subject
    rec_fnames = dir([rec_dir '/' subject_folders(s).name '/units/']);
    rec_fnames = rec_fnames(~ismember({rec_fnames(:).name},{'.','..'}));
    
    % now loop over each recording
    for r = 1:numel(rec_fnames)
        
        fprintf(['\nprocessing ' rec_fnames(r).name '\n']);
        
        % load this recording file (loads a struct called 'sdata')
        load([rec_fnames(r).folder '/' rec_fnames(r).name]);
        
        t_mids = sdata.t_mids+1;
        nrmlFRs = sdata.nrml_FRs;
        firing_rates = sdata.firing_rates;
        u_names = sdata.unit_names;
        
        % extract the behavior
        bhv = array2table(sdata.bhv);
        bhv.Properties.VariableNames = sdata.bhv_varnames;
        
        % create a grouping variable for computing trial/condition-averaged data
        state_groups = zeros(size(bhv.state));
        state_groups(bhv.state == 1 & bhv.state_type == 1) = 1;
        state_groups(bhv.state == 1 & bhv.state_type == -1) = 2;
        state_groups(bhv.state == -1 & bhv.state_type == 1) = 3;
        state_groups(bhv.state == -1 & bhv.state_type == -1) = 4;
        state1_ix = bhv.state == 1;
        state2_ix = bhv.state == -1;
        
        cue_FRs = squeeze(mean(nrmlFRs(:,51:71,:),2));

        % time points for Don
        if contains(subject_folders(s).name, 'Don')
            choice_FRs = squeeze(mean(nrmlFRs(:,87:95,:),2));

        else % time points for King
            choice_FRs = squeeze(mean(nrmlFRs(:,83:99,:),2));
        end
        
        c_FR2 = squeeze(mean(nrmlFRs(:,83:99,:),2));
        
        % initialize a table to save the selectivity info
        s_table = table;
        % loop over the neurons
        for u = 1:size(cue_FRs,2)
            u_ctr = u_ctr + 1;
            
            % do the selectivity analysis
            s_table(u,:) = CFV_assess_epochs(cue_FRs(:,u), choice_FRs(:,u), c_FR2(:,u), bhv, sdata.unit_names(u),...
                                             sdata.OFC_ix(u), {subject_folders(s).name});
            
            % get the trial-averaged data for each condition
            unit_data.subject{u_ctr} = subject_folders(s).name;
            unit_data.u_name{u_ctr} = sdata.unit_names(u);
            unit_data.OFC_ix{u_ctr} = sdata.OFC_ix(u);
            
            % get mean state firing rates
           [unit_data.state_mean_FRs{u_ctr}, gnames] = grpstats(firing_rates(:,:,u), state_groups, {'mean', 'gname'});
            
           [unit_data.state1_val_FRs{u_ctr}, s1_vals] = grpstats(firing_rates(state1_ix,:,u),...
                                                       bhv.chosenval(state1_ix), {'mean', 'gname'}); 
                                                   
           [unit_data.state2_val_FRs{u_ctr}, s2_vals] = grpstats(firing_rates(state2_ix,:,u),...
                                                       bhv.chosenval(state2_ix), {'mean', 'gname'}); 
             
        end % of looping over units in this recording
        
        s_table.cue_phase_locked = sdata.cue_phase_locking(:,1) < .05;
        s_table.choice_phase_locked = sdata.choice_phase_locking(:,1) < .05;

        % save the selectivity information into sdata
        sdata.s_table = s_table;
        
        % re-save sdata
        save([rec_fnames(r).folder '/' rec_fnames(r).name],'sdata','-v7.3');
        
        % accumulate the selectivity info
        selectivity_table = [selectivity_table ; s_table];
        
    end % of looping over each recording
    
end % of looping over subjects

fprintf('\n Saving aggregate output \n ');
unit_data.selectivity_table = selectivity_table;
unit_data.t_mids = t_mids;
save([rec_dir '/unit_data.mat'],'unit_data');

fprintf('\n all done :] \n ');