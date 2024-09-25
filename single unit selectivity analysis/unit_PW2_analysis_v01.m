% unit_selectivity_analysis.m

%-----------------------------
% This function does a sliding-window analysis of partial eta squared
%-----------------------------
rec_dir = 'C:/Users/Thomas Elston/Documents/MATLAB/Projects/CuedFlexVal/preprocessed data/';

% get subject-level directories
subject_folders = dir(rec_dir);
subject_folders = subject_folders(~ismember({subject_folders(:).name},{'.','..'}));
subject_folders = subject_folders(cell2mat({subject_folders(:).isdir}));

fprintf('\nrunning unit selectivity analysis \n');

PW2 = [];
free_PW2 = [];
forced_PW2 = [];
p_vals = [];
choice_pvals = [];
brain_area = [];
subj = [];
u_ctr = 0;
s_tbl = table;

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
        
        t_mids = (sdata.t_mids+1)/1000;
        nrmlFRs = sdata.nrml_FRs;
        firing_rates = sdata.firing_rates;
        u_names = sdata.unit_names;

        brain_area = [brain_area ; sdata.OFC_ix];
        subj = [subj; ones(size(sdata.OFC_ix))*s];

        s_tbl = [s_tbl ; sdata.s_table];
        
        % extract the behavior
        bhv = array2table(sdata.bhv);
        bhv.Properties.VariableNames = sdata.bhv_varnames;

        % drop incorrect trials
        drop_ix =  (bhv.pickedbest == 0);

        bhv(drop_ix, :) = [];
        nrmlFRs(drop_ix, :, :) = [];
        free_ix = bhv.forcedchoice == 0; 
        forced_ix = bhv.forcedchoice == 1; 
        choice_FRs = squeeze(nanmean(nrmlFRs(:, 79:95, :), 2));
        
        % create a grouping variable for computing trial/condition-averaged data
        state_ix = bhv.state;
        val_ix = bhv.chosenval;
        val_ix(bhv.chosenval == 1) = -2;
        val_ix(bhv.chosenval == 2) = -1;
        val_ix(bhv.chosenval == 3) =  1;
        val_ix(bhv.chosenval == 4) =  2;
        factor_table = table; 
        factor_table.state = state_ix; 
        factor_table.value = val_ix;

        factors = {state_ix, val_ix};
        free_factors = {state_ix(free_ix), val_ix(free_ix)};
        forced_factors = {state_ix(forced_ix), val_ix(forced_ix)};
        n_factors = numel(factors)+1;
        factor_labels = {'state', 'val', 'interaction'};

        r_PW2 = NaN(n_factors, numel(t_mids), size(nrmlFRs,3));
        r_free_PW2 = NaN(n_factors, numel(t_mids), size(nrmlFRs,3));
        r_forced_PW2 = NaN(n_factors, numel(t_mids), size(nrmlFRs,3));
        r_pvals = NaN(n_factors, numel(t_mids), size(nrmlFRs,3));

        r_choice_pvals = NaN(size(nrmlFRs,3), n_factors);


        n_times = numel(t_mids);

        % loop over the neurons
        parfor u = 1:size(nrmlFRs,3)

            [r_choice_pvals(u, :), ~] = compute_regression(choice_FRs(:, u), factor_table)

            % loop over times
            for t = 1:n_times
                [r_pvals(:, t, u), r_PW2(:, t, u)] = compute_PW2(nrmlFRs(:, t, u), factors, factor_labels);
                [~, r_free_PW2(:, t, u)] = compute_PW2(nrmlFRs(free_ix, t, u), free_factors, factor_labels);
                [~, r_forced_PW2(:, t, u)] = compute_PW2(nrmlFRs(forced_ix, t, u), forced_factors, factor_labels);


            end % of looping over times

        end % of looping over units in this recording

        PW2 = cat(3, PW2, r_PW2);
        free_PW2 = cat(3, free_PW2, r_free_PW2);
        forced_PW2 = cat(3, forced_PW2, r_forced_PW2);
        p_vals = cat(3, p_vals, r_pvals);
        choice_pvals = [choice_pvals; r_choice_pvals];

    end % of looping over each recording
    
end % of looping over subjects


k_ix = subj == 2;
d_ix = subj == 1; 
ofc_ix = brain_area == 1; 
hpc_ix = brain_area == 0;
p_vals2 = p_vals < .05;

state_p = squeeze(p_vals2(1, :,:))';
val_p = squeeze(p_vals2(2, :,:))';
stateval_p = squeeze(p_vals2(3, :,:))';

sel_ix = (sum(state_p(:,83:99), 2) >10) | (sum(state_p(:,51:67), 2) >10) |...
         (sum(stateval_p(:,83:99), 2) >10) | (sum(stateval_p(:,51:67), 2) >10);



% compute mean and CIs for each factor and brain area
[K_OFC_state_mean, K_OFC_state_ci] = GetMeanCI(squeeze(PW2(3,:, ofc_ix  & k_ix & sel_ix))'*100, 'sem');
[K_HPC_state_mean, K_HPC_state_ci] = GetMeanCI(squeeze(PW2(1,:, hpc_ix  & k_ix & sel_ix))'*100, 'sem');

[D_OFC_state_mean, D_OFC_state_ci] = GetMeanCI(squeeze(PW2(3,:, ofc_ix & d_ix & sel_ix))'*100, 'sem');
[D_HPC_state_mean, D_HPC_state_ci] = GetMeanCI(squeeze(PW2(1,:, hpc_ix & d_ix & sel_ix))'*100, 'sem');


% state props
[D_OFC_state_n, ~] = GetMeanCI(squeeze(p_vals2(3,:, ofc_ix & d_ix))'*100, 'sem');
[D_HPC_state_n, ~] = GetMeanCI(squeeze(p_vals2(1,:, hpc_ix & d_ix))'*100, 'sem');

[K_OFC_state_n, ~] = GetMeanCI(squeeze(p_vals2(3,:, ofc_ix & k_ix))'*100, 'sem');
[K_HPC_state_n, ~] = GetMeanCI(squeeze(p_vals2(1,:, hpc_ix & k_ix))'*100, 'sem');


% split by free and forced choice
[K_OFC_free_state_mean, K_OFC_free_state_ci] = GetMeanCI(squeeze(free_PW2(3,:, ofc_ix  & k_ix & sel_ix))'*100, 'sem');
[K_HPC_free_state_mean, K_HPC_free_state_ci] = GetMeanCI(squeeze(free_PW2(1,:, hpc_ix  & k_ix & sel_ix))'*100, 'sem');

[K_OFC_forced_state_mean, K_OFC_forced_state_ci] = GetMeanCI(squeeze(forced_PW2(3,:, ofc_ix  & k_ix & sel_ix))'*100, 'sem');
[K_HPC_forced_state_mean, K_HPC_forced_state_ci] = GetMeanCI(squeeze(forced_PW2(1,:, hpc_ix  & k_ix & sel_ix))'*100, 'sem');

[D_OFC_free_state_mean, D_OFC_free_state_ci] = GetMeanCI(squeeze(free_PW2(3,:, ofc_ix  & d_ix & sel_ix))'*100, 'sem');
[D_HPC_free_state_mean, D_HPC_free_state_ci] = GetMeanCI(squeeze(free_PW2(1,:, hpc_ix  & d_ix & sel_ix))'*100, 'sem');

[D_OFC_forced_state_mean, D_OFC_forced_state_ci] = GetMeanCI(squeeze(forced_PW2(3,:, ofc_ix  & d_ix & sel_ix))'*100, 'sem');
[D_HPC_forced_state_mean, D_HPC_forced_state_ci] = GetMeanCI(squeeze(forced_PW2(1,:, hpc_ix  & d_ix & sel_ix))'*100, 'sem');


n_k_OFC = squeeze(p_vals(3,:, ofc_ix & k_ix))';
n_k_HPC = squeeze(p_vals(1,:, hpc_ix & k_ix))';

n_d_OFC = squeeze(p_vals(3,:, ofc_ix & d_ix))';
n_d_HPC = squeeze(p_vals(1,:, hpc_ix & d_ix))';

sum_k_OFC = sum(sum(n_k_OFC(:,83:99), 2) >4);
sum_k_HPC = sum(sum(n_k_HPC(:,83:99), 2) >4);

sum_d_OFC = sum(sum(n_d_OFC(:,83:99), 2) >4);
sum_d_HPC = sum(sum(n_d_HPC(:,83:99), 2) >4);

[tbl,chi2stat,pval] = chiSquareWithFrequencies_v02(sum_d_OFC, sum(ofc_ix & d_ix), ...
                                                   sum_d_HPC, sum(hpc_ix & d_ix))


cmap = cbrewer('qual', 'Accent', 8);
hpc_c = cmap(8,:);
ofc_c = cmap(5,:);

figure;
subplot(1,2,1)
hold on
plot(t_mids, smoothdata(K_OFC_state_n, "movmean", 5), 'color', ofc_c, 'LineWidth', 2);
plot(t_mids, smoothdata(K_HPC_state_n, "movmean", 5), 'color', hpc_c, 'LineWidth', 2);
xlim([-1, 1])
title('Subject K')
ylabel('% State-Selective Units')
axis square


subplot(1,2,2)
hold on
plot(t_mids, smoothdata(D_OFC_state_n, "movmean", 5), 'color', ofc_c, 'LineWidth', 2);
plot(t_mids, smoothdata(D_HPC_state_n, "movmean", 5), 'color', hpc_c, 'LineWidth', 2);
xlim([-1, 1])
title('Subject D')
axis square



figure;
subplot(1,2,1)
hold on
shadedErrorBar(t_mids, K_OFC_state_mean, K_OFC_state_ci, 'lineprops',{'color', ofc_c, 'LineWidth', 2});
shadedErrorBar(t_mids, K_HPC_state_mean+.2, K_HPC_state_ci, 'lineprops',{'color', hpc_c, 'LineWidth', 2});
ylim([0, 5])
xlim([-1, 1])
title('Subject K')
ylabel('Partial Omega Squared')
axis square

subplot(1, 2, 2)
hold on
shadedErrorBar(t_mids, D_OFC_state_mean, D_OFC_state_ci, 'lineprops',{'color', ofc_c, 'LineWidth', 2});
shadedErrorBar(t_mids, D_HPC_state_mean, D_HPC_state_ci, 'lineprops',{'color', hpc_c, 'LineWidth', 2});
xlim([-1, 1])
ylim([0, 4])
axis square
title('Subject D')


figure;
subplot(2,2,1)
hold on
shadedErrorBar(t_mids, K_OFC_free_state_mean, K_OFC_free_state_ci, 'lineprops',{'color', ofc_c, 'LineWidth', 2});
shadedErrorBar(t_mids, K_HPC_free_state_mean, K_HPC_free_state_ci, 'lineprops',{'color', hpc_c, 'LineWidth', 2});
ylim([0, 6])
xlim([-1, 1])
title('Subject K')
axis square

subplot(2,2,2)
hold on
shadedErrorBar(t_mids, D_OFC_free_state_mean, D_OFC_free_state_ci, 'lineprops',{'color', ofc_c, 'LineWidth', 2});
shadedErrorBar(t_mids, D_HPC_free_state_mean, D_HPC_free_state_ci, 'lineprops',{'color', hpc_c, 'LineWidth', 2});
ylim([0, 5])
xlim([-1, 1])
title('Subject D')
axis square

subplot(2,2,3)
hold on
shadedErrorBar(t_mids, K_OFC_forced_state_mean, K_OFC_forced_state_ci, 'lineprops',{'color', ofc_c, 'LineWidth', 2});
shadedErrorBar(t_mids, K_HPC_forced_state_mean, K_HPC_forced_state_ci, 'lineprops',{'color', hpc_c, 'LineWidth', 2});
ylim([0, 7])
xlim([-1, 1])
axis square

subplot(2,2,4)
hold on
shadedErrorBar(t_mids, D_OFC_forced_state_mean, D_OFC_forced_state_ci, 'lineprops',{'color', ofc_c, 'LineWidth', 2});
shadedErrorBar(t_mids, D_HPC_forced_state_mean, D_HPC_forced_state_ci, 'lineprops',{'color', hpc_c, 'LineWidth', 2});
ylim([0, 7])
xlim([-1, 1])
axis square



%----------------------
d_ofc = ofc_ix & d_ix; 
n_d_ofc = sum(d_ofc);
d_hpc = hpc_ix & d_ix; 
n_d_hpc = sum(d_hpc);
d_OFC_state_only = sum(choice_pvals(d_ofc, 1) < .05 & choice_pvals(d_ofc, 2) > .05 & choice_pvals(d_ofc, 3) > .05);
d_OFC_val_only =   sum(choice_pvals(d_ofc, 1) > .05 & choice_pvals(d_ofc, 2) < .01 & choice_pvals(d_ofc, 3) > .05);
d_OFC_state_val =  sum(                               choice_pvals(d_ofc, 2) > .01 & choice_pvals(d_ofc, 3) < .05);

d_HPC_state_only = sum(choice_pvals(d_hpc, 1) < .05 & choice_pvals(d_hpc, 2) > .05 & choice_pvals(d_hpc, 3) > .05);
d_HPC_val_only =   sum(choice_pvals(d_hpc, 1) > .05 & choice_pvals(d_hpc, 2) < .05 & choice_pvals(d_hpc, 3) > .01);
d_HPC_state_val =  sum(choice_pvals(d_hpc, 1) > .05 & choice_pvals(d_hpc, 2) > .05 & choice_pvals(d_hpc, 3) < .05);


don_props = [d_HPC_state_only/n_d_hpc, d_HPC_val_only/n_d_hpc, d_HPC_state_val/n_d_hpc;
             d_OFC_state_only/n_d_ofc, d_OFC_val_only/n_d_ofc, d_OFC_state_val/n_d_ofc];


k_ofc = ofc_ix & k_ix; 
n_k_ofc = sum(k_ofc);
k_hpc = hpc_ix & k_ix; 
n_k_hpc = sum(k_hpc);
k_OFC_state_only = sum(choice_pvals(k_ofc, 1) < .05 & choice_pvals(k_ofc, 2) > .05 & choice_pvals(k_ofc, 3) > .05);
k_OFC_val_only =   sum(choice_pvals(k_ofc, 1) > .05 & choice_pvals(k_ofc, 2) < .01 & choice_pvals(k_ofc, 3) > .05);
k_OFC_state_val =  sum(                               choice_pvals(k_ofc, 2) > .01 & choice_pvals(k_ofc, 3) < .05);

k_HPC_state_only = sum(choice_pvals(k_hpc, 1) < .05 & choice_pvals(k_hpc, 2) > .05 & choice_pvals(k_hpc, 3) > .05);
k_HPC_val_only =   sum(choice_pvals(k_hpc, 1) > .05 & choice_pvals(k_hpc, 2) < .05 & choice_pvals(k_hpc, 3) > .01);
k_HPC_state_val =  sum(choice_pvals(k_hpc, 1) > .05 & choice_pvals(k_hpc, 2) > .05 & choice_pvals(k_hpc, 3) < .05);



king_props = [k_HPC_state_only/n_k_hpc, k_HPC_val_only/n_k_hpc, k_HPC_state_val/n_k_hpc;
              k_OFC_state_only/n_k_ofc, k_OFC_val_only/n_k_ofc, k_OFC_state_val/n_k_ofc];

figure; 
subplot(1,2,1)
bar(king_props' * 100)
xlim([.6, 3.4])
ylim([0, 21])
axis square
title('Subject K')

subplot(1,2,2)
bar(don_props' * 100)
xlim([.6, 3.4])
axis square
title('Subject D')



% comparing pure value across brain areas
[~,d_val_X2,d_val_pval] = chiSquareWithFrequencies_v02(d_HPC_val_only, n_d_hpc, ...
                                                       d_OFC_val_only, n_d_ofc);

% comparing state value interaction across brain areas
[~,d_stateval_X2,d_stateval_pval] = chiSquareWithFrequencies_v02(d_HPC_state_val, n_d_hpc, ...
                                                                 d_OFC_state_val, n_d_ofc);
% comparing pure value against state-value within brain areas
[~,d_hpc_val_X2,d_hpc_val_pval] = chiSquareWithFrequencies_v02(d_HPC_val_only, n_d_hpc, ...
                                                               d_HPC_state_val, n_d_hpc);

[~,d_ofc_val_X2,d_ofc_val_pval] = chiSquareWithFrequencies_v02(d_OFC_val_only, n_d_ofc, ...
                                                               d_OFC_state_val, n_d_ofc);











% compare pure value across brain areas
[~,k_val_X2,k_val_pval] = chiSquareWithFrequencies_v02(k_HPC_val_only, n_k_hpc, ...
                                                       k_OFC_val_only, n_k_ofc);

% compare state-value interaction across brain areas
[~,k_stateval_X2,k_stateval_pval] = chiSquareWithFrequencies_v02(k_HPC_state_val, n_k_hpc, ...
                                                                 k_OFC_state_val, n_k_ofc);

% comparing pure value against state-value within brain areas
[~,k_hpc_val_X2,k_hpc_val_pval] = chiSquareWithFrequencies_v02(k_HPC_val_only, n_k_hpc, ...
                                                               k_HPC_state_val, n_k_hpc);

[~,k_ofc_val_X2,k_ofc_val_pval] = chiSquareWithFrequencies_v02(k_OFC_val_only, n_k_ofc, ...
                                                               k_OFC_state_val, n_k_ofc);












fprintf('\n all done :] \n ');