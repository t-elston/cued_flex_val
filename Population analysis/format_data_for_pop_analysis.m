%% format_data_for_pop_analysis.m

% High level script for puttind data into an easy-to-use format for analyzing psuedo-populations. 

% Thomas Elston
% 8 November 2023

%--------------------------------------------------------------
% Step 1: format data into "super arrays"
% Step 2: create psuedopopulations and create a bunch of bootstraps
%
% *** This function needs to be run uniquely for each subject
%--------------------------------------------------------------

% define some paths for loading/saving
% where are the units for this subject?
datadir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\preprocessed data\Don\units/';

% where do we want to save the "super arrays" for this subject?
super_array_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\preprocessed data\Don\population\super_arrays/';

% where do we want to save the bootstraps for this subject?
bootdir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\preprocessed data\Don\population\boot_PCs/';



%% STEP 1 - organize the firing rates from individual single units into one huge array of size n_units x n_time_steps.
% It also generates a corresponding table with the relevant features of each trial.

rec_files = dir([datadir '*.mat']);

% subselect and save the data from each brain area
OFC_data=struct;
OFC_data.zFRs=[];
OFC_data.u_names=[];
OFC_data.bhv=[];
OFC_data.n_trials_per_choice_cond=[];
OFC_data.brain_area = 'OFC';


HPC_data=struct;
HPC_data.zFRs=[];
HPC_data.u_names=[];
HPC_data.bhv=[];
HPC_data.n_trials_per_choice_cond=[];
HPC_data.brain_area = 'HPC';


for f = 1:numel(rec_files)
    
    mat_fname = [datadir rec_files(f).name];
    
    fprintf('\n')  
    disp(['***loading ' rec_files(f).name]) 
    load(mat_fname);
    
    % make table of the behavior and count number of trials per condition
    bhv2 = array2table(sdata.bhv,'VariableNames',sdata.bhv_varnames);

    firing_rates = sdata.nrml_FRs;
    t_mids = sdata.t_mids;
    
    [n_trials, n_times, n_units] = size(firing_rates);
    
    % now start constructing the arrays for each unit in the file
    for u = 1:n_units
        
        % was this a HPC or OFC unit?
        if sdata.OFC_ix(u) == 1
            
            OFC_data.zFRs = [OFC_data.zFRs ; firing_rates(:,:,u)];
            OFC_data.u_names = [OFC_data.u_names ; repmat(sdata.unit_names(u),n_trials,1)]; 
            OFC_data.bhv = [OFC_data.bhv ; sdata.bhv];

            % ensure there are enough trials/condition during the choice for this neuron
            OFC_data.n_trials_per_choice_cond= cat(2, OFC_data.n_trials_per_choice_cond,... 
                              [ sum(bhv2.state == 1 & bhv2.chosenval == 1) ; 
                                sum(bhv2.state == 1 & bhv2.chosenval == 2) ;
                                sum(bhv2.state == 1 & bhv2.chosenval == 3) ;
                                sum(bhv2.state == 1 & bhv2.chosenval == 4) ;
                                sum(bhv2.state == -1 & bhv2.chosenval == 1) ; 
                                sum(bhv2.state == -1 & bhv2.chosenval == 2) ;
                                sum(bhv2.state == -1 & bhv2.chosenval == 3) ;
                                sum(bhv2.state == -1 & bhv2.chosenval == 4)]) ;
           
        else
            
            HPC_data.zFRs = [HPC_data.zFRs ; firing_rates(:,:,u)];
            HPC_data.u_names = [HPC_data.u_names ; repmat(sdata.unit_names(u),n_trials,1)]; 
            HPC_data.bhv = [HPC_data.bhv ; sdata.bhv];
            
            % ensure there are enough trials/condition during the choice for this neuron
            HPC_data.n_trials_per_choice_cond= cat(2, HPC_data.n_trials_per_choice_cond,... 
                              [ sum(bhv2.state == 1 & bhv2.chosenval == 1) ; 
                                sum(bhv2.state == 1 & bhv2.chosenval == 2) ;
                                sum(bhv2.state == 1 & bhv2.chosenval == 3) ;
                                sum(bhv2.state == 1 & bhv2.chosenval == 4) ;
                                sum(bhv2.state == -1 & bhv2.chosenval == 1) ; 
                                sum(bhv2.state == -1 & bhv2.chosenval == 2) ;
                                sum(bhv2.state == -1 & bhv2.chosenval == 3) ;
                                sum(bhv2.state == -1 & bhv2.chosenval == 4)]) ;
                            
        end % of determining whether this is a HPC or OFC unit
               
    end % of looping over units
    
end % of looping over files

OFC_data.bhv_varnames = sdata.bhv_varnames;
OFC_data.ts = t_mids;

HPC_data.bhv_varnames = sdata.bhv_varnames;
HPC_data.ts = t_mids;

% now let's save the data
save([super_array_dir, 'OFC_popdata.mat'],'-struct','OFC_data','-v7.3');
save([super_array_dir, 'HPC_popdata.mat'],'-struct','HPC_data','-v7.3');


fprintf('\n')  
disp('Finished making super arrays.');

% now clear the variables and let's move to step 2 to create the bootstraps
clearvars -except datadir super_array_dir bootdir
clc

%% STEP 2 - making bootstraps


f_data = dir([super_array_dir '*.mat']);

n_boots = 1000;

for f = 1:numel(f_data)
    
    % load data
    load([super_array_dir f_data(f).name]);
    
    bhv = array2table(bhv,'VariableNames',bhv_varnames);
    u_ids = unique(u_names);
    n_units = numel(u_ids);
    n_times = size(zFRs,2);

    % find relevant epochs (in ms, relative to when choice options appear)
    [~, cue_on] = min(abs(ts - -700));
    [~, cue_off] = min(abs(ts - -200));
    [~, choice_on] = min(abs(ts - 100));
    [~, choice_off] = min(abs(ts - 500));
    
    hit_ix = bhv.pickedbest == 1;
    
    % get mean firing rates during cue and choice periods
    % cue_data = nanmean(zFRs(:, cue_on : cue_off), 2);
    choice_data =nanmean(zFRs(:, choice_on : choice_off), 2);
    
    choice_STDs=[];
    cue_FRs=[];
    choice_FRs=[];
    
    % define single PC spaces for cue and choice
    for u = 1:n_units
        
        u_ix = contains(u_names, u_ids(u));
                           
        % cue_FRs(:, u) = zscore([nanmean(cue_data(u_ix & bhv.state == 1)) ;
        %                         nanmean(cue_data(u_ix & bhv.state == -1)) ;
        %                         nanmean(cue_data(u_ix & bhv.state_type == 1)) ;
        %                         nanmean(cue_data(u_ix & bhv.state_type == -1))]);                    

       choice_FRs(:, u) = zscore([ ;
            nanmean(choice_data(u_ix & bhv.chosenval == 1 & bhv.state == 1));
            nanmean(choice_data(u_ix & bhv.chosenval == 2 & bhv.state == 1));
            nanmean(choice_data(u_ix & bhv.chosenval == 3 & bhv.state == 1));
            nanmean(choice_data(u_ix & bhv.chosenval == 4 & bhv.state == 1));
            nanmean(choice_data(u_ix & bhv.chosenval == 1 & bhv.state == -1));
            nanmean(choice_data(u_ix & bhv.chosenval == 2 & bhv.state == -1));
            nanmean(choice_data(u_ix & bhv.chosenval == 3 & bhv.state == -1));
            nanmean(choice_data(u_ix & bhv.chosenval == 4 & bhv.state == -1));
            ]);
         
    end % of looping over neurons
    

    % compute a common PC space to project all bootstraps into
    [choice_coeff, choice_score, choice_latent, ~, choice_expvar, choice_mu] = pca(choice_FRs);
    n_choice_PCs = 5;

    % [cue_coeff, cue_score, cue_latent, ~, cue_expvar, cue_mu] = pca(cue_FRs);
    % n_cue_PCs = 3;
       
    cue_PCs=[];
    choice_PCs=[];
    
    pw = PoolWaitbar(n_boots, 'bootstrapping...');
    parfor b = 1:n_boots
        increment(pw);
        
        % b_cue_FRs = [];
        b_choice_FRs = [];
        %fprintf(['\n region: ' f_data(f).name(1:3) ', b = ' num2str(b) ' / ' num2str(n_boots)]);
        
        for u = 1:n_units
            
            u_ix = contains(u_names, u_ids(u));
            
            % % grab 80% of trials for each condition during cue
            % c1_t1 = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.state_type == 1,.8);
            % c1_t2 = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.state_type == -1,.8);
            % c2_t1 = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.state_type == 1,.8);
            % c2_t2 = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.state_type == -1,.8);
            
            % grab 80% of the trials for each condition during choice
            c1_v1 = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.chosenval == 1,.8);
            c1_v2 = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.chosenval == 2,.8);
            c1_v3 = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.chosenval == 3,.8);
            c1_v4 = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.chosenval == 4,.8);
            c2_v1 = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.chosenval == 1,.8);
            c2_v2 = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.chosenval == 2,.8);
            c2_v3 = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.chosenval == 3,.8);
            c2_v4 = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.chosenval == 4,.8);
                        
            % b_cue_FRs(:, u) = zscore([nanmean(cue_data(c1_t1)) ;
            %     nanmean(cue_data(c1_t2)) ;
            %     nanmean(cue_data(c2_t1)) ;
            %     nanmean(cue_data(c2_t2))]);
            
            b_choice_FRs(:, u) = zscore([nanmean(choice_data(c1_v1)) ;
                nanmean(choice_data(c1_v2)) ;
                nanmean(choice_data(c1_v3)) ;
                nanmean(choice_data(c1_v4)) ;
                nanmean(choice_data(c2_v1)) ;
                nanmean(choice_data(c2_v2)) ;
                nanmean(choice_data(c2_v3)) ;
                nanmean(choice_data(c2_v4))]);
        
        end % of looping over neurons
        
        % project this bootstrap's firing rates into a common PC space aggregate
        % [cue_PCs(:,:,b)] = ProjectIntoPCspace_v02(b_cue_FRs,cue_mu,cue_coeff,n_cue_PCs);
        [choice_PCs(:,:,b)] = ProjectIntoPCspace_v02(b_choice_FRs,choice_mu,choice_coeff,n_choice_PCs);
        
    end % of looping over bootstraps
    delete(pw);
    
    if contains(brain_area, 'HPC')
        % HPC_cue_PCs = cue_PCs;
        HPC_choice_PCs = choice_PCs;
    else
        % OFC_cue_PCs = cue_PCs;
        OFC_choice_PCs = choice_PCs;
    end


    
end % of looping over files names

cue_ctx_labels = [1;1;2;2];
cue_type_labels = [1;2;1;2];
choice_ctx_labels = [1;1;1;1;2;2;2;2];
choice_val_labels = [1;2;3;4;1;2;3;4];


save([bootdir 'popdata'], 'HPC_choice_PCs', 'OFC_choice_PCs', 'choice_val_labels', 'choice_ctx_labels')






















