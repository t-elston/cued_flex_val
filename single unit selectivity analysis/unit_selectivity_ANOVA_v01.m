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
aov_table = table;
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
        
        % create factors for ANOVA analysis
        factors = {bhv.state, bhv.chosenval};
        fix_
        cue_FRs = squeeze(mean(nrmlFRs(:,51:71,:),2));
        choice_FRs = squeeze(mean(nrmlFRs(:,83:99,:),2));

        
        % initialize a table to save the selectivity info
        s_table = table;
        % loop over the neurons
        for u = 1:size(cue_FRs,2)
            u_ctr = u_ctr + 1;
           
             
        end % of looping over units in this recording
      
        
    end % of looping over each recording
    
end % of looping over subjects

fprintf('\n Saving aggregate output \n ');
unit_data.selectivity_table = selectivity_table;
unit_data.t_mids = t_mids;
save([rec_dir '/unit_data.mat'],'unit_data');

fprintf('\n all done :] \n ');