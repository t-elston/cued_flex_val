% aaa_CFV_pop_trajectory_main_v02.m

datadir = 'C:\Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/population_data/';

f_data = dir([datadir '*.mat']);

n_boots = 1000;

for f = 1:numel(f_data)
    
    % load data
    load([datadir f_data(f).name]);
    
    bhv = array2table(bhv,'VariableNames',bhv_varnames);
    u_ids = unique(u_names);
    n_units = numel(u_ids);
    n_times = size(zFRs,2);
    cue_on = find(ts == -700);
    cue_off = find(ts == -200);
    choice_on = find(ts == 0);
    choice_off = find(ts == 500);
    
    hit_ix = bhv.pickedbest == 1;
    
    % get mean firing rates during cue and choice periods
    cue_data = mean(zFRs(:, cue_on : cue_off), 2);
    choice_data = mean(zFRs(:, choice_on : choice_off), 2);
    
    choice_STDs=[];
    cue_FRs=[];
    choice_FRs=[];
    
    % define single PC spaces for cue and choice
    for u = 1:n_units
        
        u_ix = contains(u_names, u_ids(u));
                           
        cue_FRs(:, u) = zscore([nanmean(cue_data(u_ix & bhv.state == 1)) ;
                                nanmean(cue_data(u_ix & bhv.state == -1)) ;
                                nanmean(cue_data(u_ix & bhv.state_type == 1)) ;
                                nanmean(cue_data(u_ix & bhv.state_type == -1))]);                    

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
%         
    end % of looping over neurons
    

    % compute and store the PCs
    [cue_coeff, cue_score, cue_latent, ~, cue_expvar, cue_mu] = pca(cue_FRs);
    [choice_coeff, choice_score, choice_latent, ~, choice_expvar, choice_mu] = pca(choice_FRs); 
    n_cue_PCs = 3;
    n_choice_PCs = 5;
    
      % do ICA
     cue_mdl = rica(cue_FRs, 2);
%      choice_mdl = rica(choice_FRs, 3);
%     
    cue_PCs=[];
    choice_PCs=[];
    
    pw = PoolWaitbar(n_boots, 'bootstrapping...');
    parfor b = 1:n_boots
        increment(pw);
        
        b_cue_FRs = [];
        b_choice_FRs = [];
        %fprintf(['\n region: ' f_data(f).name(1:3) ', b = ' num2str(b) ' / ' num2str(n_boots)]);
        
        for u = 1:n_units
            
            u_ix = contains(u_names, u_ids(u));
            
            % grab 80% of trials for each condition during cue
            c1_t1 = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.state_type == 1,.8);
            c1_t2 = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.state_type == -1,.8);
            c2_t1 = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.state_type == 1,.8);
            c2_t2 = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.state_type == -1,.8);
            
            % grab 80% of the trials for each condition during choice
            c1_v1 = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.chosenval == 1,.8);
            c1_v2 = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.chosenval == 2,.8);
            c1_v3 = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.chosenval == 3,.8);
            c1_v4 = get_random_proportion_v02(u_ix & bhv.state == 1 & bhv.chosenval == 4,.8);
            c2_v1 = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.chosenval == 1,.8);
            c2_v2 = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.chosenval == 2,.8);
            c2_v3 = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.chosenval == 3,.8);
            c2_v4 = get_random_proportion_v02(u_ix & bhv.state == -1 & bhv.chosenval == 4,.8);
                        
            b_cue_FRs(:, u) = zscore([nanmean(cue_data(c1_t1)) ;
                nanmean(cue_data(c1_t2)) ;
                nanmean(cue_data(c2_t1)) ;
                nanmean(cue_data(c2_t2))]);
            
            b_choice_FRs(:, u) = zscore([nanmean(choice_data(c1_v1)) ;
                nanmean(choice_data(c1_v2)) ;
                nanmean(choice_data(c1_v3)) ;
                nanmean(choice_data(c1_v4)) ;
                nanmean(choice_data(c2_v1)) ;
                nanmean(choice_data(c2_v2)) ;
                nanmean(choice_data(c2_v3)) ;
                nanmean(choice_data(c2_v4))]);
        
        end % of looping over neurons
        
        % project this bootstrap's PCs and aggregate
%          cue_PCs(:,:,b) = transform(cue_mdl, b_cue_FRs);
%          choice_PCs(:,:,b) = transform(choice_mdl, b_choice_FRs);

        [cue_PCs(:,:,b)] = ProjectIntoPCspace_v02(b_cue_FRs,cue_mu,cue_coeff,n_cue_PCs);
        [choice_PCs(:,:,b)] = ProjectIntoPCspace_v02(b_choice_FRs,choice_mu,choice_coeff,n_choice_PCs);
        
    end % of looping over bootstraps
    delete(pw);
    
    if f == 1
        HPC_cue_PCs = cue_PCs;
        HPC_choice_PCs = choice_PCs;
    else
        OFC_cue_PCs = cue_PCs;
        OFC_choice_PCs = choice_PCs;
    end
    
end % of looping over files names


cue_ctx_labels = [1;1;2;2];
cue_type_labels = [1;2;1;2];
choice_ctx_labels = [1;1;1;1;2;2;2;2];
choice_val_labels = [1;2;3;4;1;2;3;4];




% manually save the data and do analysis in another function



