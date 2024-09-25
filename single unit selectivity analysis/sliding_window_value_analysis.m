% sliding_window_value_analysis.m

%-----------------------------
% This function does a sliding-window analysis of value coding
%-----------------------------
rec_dir = 'C:/Users/Thomas Elston/Documents/MATLAB/Projects/CuedFlexVal/preprocessed data/';

% get subject-level directories
subject_folders = dir(rec_dir);
subject_folders = subject_folders(~ismember({subject_folders(:).name},{'.','..'}));
subject_folders = subject_folders(cell2mat({subject_folders(:).isdir}));

fprintf('\nrunning unit selectivity analysis \n');

a_p_vals = [];
b_p_vals = [];
a_betas = [];
b_betas = [];
choice_betas = [];
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
        n_times = numel(t_mids);
        n_units = numel(u_names);

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
        
        % create value regressors
        state_a_val = bhv.chosenval(bhv.state == -1);
        state_b_val = bhv.chosenval(bhv.state == 1);
        a_regressor = [ones(size(state_a_val)), state_a_val];
        b_regressor = [ones(size(state_b_val)), state_b_val];

        % extract firing rates for the states
        state_a_frs = nrmlFRs(bhv.state == -1, :,:);
        state_b_frs = nrmlFRs(bhv.state == 1, :,:);

        state_a_choice_fr = squeeze(nanmean(state_a_frs(:,83:95,:), 2));
        state_b_choice_fr = squeeze(nanmean(state_b_frs(:,83:95,:), 2));

        % initialize arrays to accumulate results into
        r_a_p_vals = NaN(n_units, n_times);
        r_b_p_vals = NaN(n_units, n_times);
        r_a_betas = NaN(n_units, n_times);
        r_b_betas = NaN(n_units, n_times);

        r_choice_pvals = NaN(n_units, 2); % first col = state 1, second col = state 2
        r_choice_betas = NaN(n_units, 2); % first col = state 1, second col = state 2


        % loop over the neurons
        parfor u = 1:size(nrmlFRs,3)

            [c_a_betas, ~, ~, ~, c_a_stats] = regress(state_a_choice_fr(:,u), a_regressor);
            [c_b_betas, ~, ~, ~, c_b_stats] = regress(state_b_choice_fr(:,u), b_regressor);

            r_choice_betas(u,:) = [c_a_betas(2), c_b_betas(2)];
            r_choice_pvals(u,:) = [c_a_stats(3), c_b_stats(3)];

            % loop over times
            for t = 1:n_times

                [t_a_betas, ~, ~, ~, a_stats] = regress(state_a_frs(:,t, u), a_regressor);
                [t_b_betas, ~, ~, ~, b_stats] = regress(state_b_frs(:,t, u), b_regressor);

                r_a_betas(u, t) = t_a_betas(2);
                r_a_p_vals(u, t) = a_stats(3);
                r_b_betas(u, t) = t_b_betas(2);
                r_b_p_vals(u, t) = b_stats(3);

            end % of looping over times

        end % of looping over units in this recording

        a_betas = [a_betas ; r_a_betas];
        a_p_vals = [a_p_vals ; r_a_p_vals];
        b_betas = [b_betas ; r_b_betas];
        b_p_vals = [b_p_vals ; r_b_p_vals];

        choice_betas = [choice_betas; r_choice_betas];
        choice_pvals = [choice_pvals; r_choice_pvals];

    end % of looping over each recording
    
end % of looping over subjects

% generate some indices
k_ix = subj == 2;
d_ix = subj == 1; 
ofc_ix = brain_area == 1; 
hpc_ix = brain_area == 0;

% now we need to regroup the data into preferred and non-preferred states
a_sig = choice_pvals(:, 1) < .05;
b_sig = choice_pvals(:, 2) < .05;
pref_state_a = (abs(choice_betas(:,1)) - abs(choice_betas(:, 2)) > 0) & a_sig;
pref_state_b = (abs(choice_betas(:,2)) - abs(choice_betas(:, 1)) > 0) & b_sig;

pref_state_p_vals = [a_p_vals(pref_state_a, :) ; b_p_vals(pref_state_b, :)] < .05;
pref_state_p_betas = abs([a_betas(pref_state_a, :) ; b_betas(pref_state_b, :)]);
pref_subj = [subj(pref_state_a) ; subj(pref_state_b)];
pref_area = [brain_area(pref_state_a) ; brain_area(pref_state_b)];

[k_hpc_beta_mean, k_hpc_beta_ci] = GetMeanCI(pref_state_p_betas(pref_subj ==2 & pref_area == 0,:), 'bootstrap');
[k_ofc_beta_mean, k_ofc_beta_ci] = GetMeanCI(pref_state_p_betas(pref_subj ==2 & pref_area == 1,:), 'bootstrap');
[d_hpc_beta_mean, d_hpc_beta_ci] = GetMeanCI(pref_state_p_betas(pref_subj ==1 & pref_area == 0,:), 'bootstrap');
[d_ofc_beta_mean, d_ofc_beta_ci] = GetMeanCI(pref_state_p_betas(pref_subj ==1 & pref_area == 1,:), 'bootstrap');

% now we need to get the latencies of value encoding
[~, choice_on] = min(abs(t_mids - .1)); 
[~, choice_off] = min(abs(t_mids - .4)); 

choice_sigs = pref_state_p_vals(:, choice_on:choice_off);

[~,latencies]=max(choice_sigs,[],2);
latency_ix = choice_on + latencies;
latency_betas = pref_state_p_betas(:, latency_ix(:));
latencies = (latencies*25) + 50; % convert to milliseconds from samples

cmap = cbrewer('qual', 'Accent', 8);
hpc_c = cmap(8,:);
ofc_c = cmap(5,:);

[p, table] = anova1(latencies(pref_subj ==2), pref_area(pref_subj ==2), 'covar', latency_betas(pref_subj ==2));
[p, table] = anova1(latencies(pref_subj ==1), pref_area(pref_subj ==1), 'covar', latency_betas(pref_subj ==1));



figure; 
subplot(2,2,1)
hold on
plot(t_mids+.05, mean(pref_state_p_vals(pref_subj ==2 & pref_area == 0,:)*100), 'color', hpc_c, 'LineWidth', 2)
plot(t_mids+.05, mean(pref_state_p_vals(pref_subj ==2 & pref_area == 1,:)*100), 'color', ofc_c, 'LineWidth', 2)
xlim([-1,1])
ylim([0, 90])
title('Subject K')
ylabel('Prop. Sig. Units')
axis square

subplot(2,2,2)
hold on
plot(t_mids+.05, mean(pref_state_p_vals(pref_subj ==1 & pref_area == 0,:)*100), 'color', hpc_c, 'LineWidth', 2)
plot(t_mids+.05, mean(pref_state_p_vals(pref_subj ==1 & pref_area == 1,:)*100), 'color', ofc_c, 'LineWidth', 2)
xlim([-1,1])
ylim([0, 90])
title('Subject D')
axis square


subplot(2,2,3)
hold on
shadedErrorBar(t_mids+.05, k_hpc_beta_mean, k_hpc_beta_ci, 'LineProps', {'color', hpc_c, 'LineWidth', 2})
shadedErrorBar(t_mids+.05, k_ofc_beta_mean, k_ofc_beta_ci, 'LineProps', {'color', ofc_c, 'LineWidth', 2})
xlim([-1,1])
ylim([0, .16])
ylabel('abs(beta of preferred state)')
axis square

subplot(2,2,4)
hold on
shadedErrorBar(t_mids+.05, d_hpc_beta_mean, d_hpc_beta_ci, 'LineProps', {'color', hpc_c, 'LineWidth', 2})
shadedErrorBar(t_mids+.05, d_ofc_beta_mean, d_ofc_beta_ci, 'LineProps', {'color', ofc_c, 'LineWidth', 2})
xlim([-1,1])
ylim([0, .3])
axis square









