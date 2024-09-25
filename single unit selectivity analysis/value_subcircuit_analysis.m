% value_subcircuit_analysis.m

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
subcircuit_tbl = table;


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
               
        choice_FRs = squeeze(mean(nrmlFRs(:,83:99,:),2));
        stateA_ix = bhv.state == 1; 
        stateB_ix = bhv.state == -1; 
        val = bhv.chosenval;
        
        % center the value code
        val(bhv.chosenval==4) = 2; 
        val(bhv.chosenval==3) = 1; 
        val(bhv.chosenval==2) = -1; 
        val(bhv.chosenval==1) = -2; 
        
        % set the predictions of the alternative state to 0
        stateA_val = val; stateA_val(stateB_ix) = 0;
        stateB_val = val; stateB_val(stateA_ix) = 0;
        
        % fill out a table to run the regressions on
        reg_tbl = table;
        reg_tbl.all_val = val;
        reg_tbl.stateA_val = stateA_val;
        reg_tbl.stateB_val = stateB_val;
        reg_tbl.state_x_val = val.*bhv.state;
        
       
        % initialize a table to save the selectivity info
        circuit_table = table;
        % loop over the neurons
        for u = 1:size(choice_FRs,2)
            
            % add the firing rates to the table
            reg_tbl.FR = choice_FRs(:,u);
            
            % run two regressions - one with all value, one with state-dependent value
            all_mdl = fitglm(reg_tbl, 'FR ~ all_val');
            stateval_mdl = fitglm(reg_tbl, 'FR ~ stateA_val + stateB_val');

            circuit_table.u_name(u) = sdata.unit_names(u);
            circuit_table.subject(u) = s;
            circuit_table.ofc(u) = sdata.OFC_ix(u);
            circuit_table.all_AIC(u) = all_mdl.ModelCriterion.AIC;
            circuit_table.stateval_AIC(u) = stateval_mdl.ModelCriterion.AIC;
            circuit_table.AIC_diff(u) = stateval_mdl.ModelCriterion.AIC - all_mdl.ModelCriterion.AIC;

        end % of looping over units in this recording
        

        % save the selectivity information into sdata
        sdata.circuit_table = circuit_table;
        
        % re-save sdata
        save([rec_fnames(r).folder '/' rec_fnames(r).name],'sdata','-v7.3');
        
        % accumulate the selectivity info
        subcircuit_tbl = [subcircuit_tbl ; circuit_table];
        
    end % of looping over each recording
    
end % of looping over subjects

fprintf('\n Saving aggregate output \n ');

save([rec_dir '/subcircuit_tbl.mat'],'subcircuit_tbl');

fprintf('\n all done :] \n ');
















