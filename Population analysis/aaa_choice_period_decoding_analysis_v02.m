% Choice period decoding analysis

base_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\preprocessed data/';

subj_names = {'Don', 'King'};
brain_areas = {'HPC', 'OFC'};


% details of an 8 way decoder
choice_trial_types = [ 1, 1, 1, 1,-1,-1,-1,-1; % states
                       1, 2, 3, 4, 1, 2, 3, 4]'; % values

% how many times should we run the decoder?
n_boots = 1000;

gcca_results = struct;

% loop over each subject
for s = 1:numel(subj_names)

    % get names of files in this directory
    f_names_struct = dir(fullfile(base_dir, subj_names{s}, 'population', 'super_arrays'));
    bytes = [f_names_struct.bytes];
    f_names = {f_names_struct(bytes > 0).name};

    % loop over each brain area
    for a = 1:numel(brain_areas)

        % find the name of the file associated with this brain area
        this_fname = f_names{contains(f_names, brain_areas{a})};

        % load the data of this brain area
        data = load([f_names_struct(1).folder '/' this_fname]);

        % extract the behavior for this brain area
        bhv = array2table(data.bhv);
        bhv.Properties.VariableNames = data.bhv_varnames;

        % get the mean firing rate over the choice period
        [~, win_start_ix] = min(abs(data.ts - 50));
        [~, win_end_ix] = min(abs(data.ts - 550));

        % mean firing rate for every trial of every neuron
        mean_fr = nanmean(data.zFRs(:, win_start_ix: win_end_ix), 2);

        % get the names of the neurons
        unit_per_trial = data.u_names;
        neuron_names = unique(unit_per_trial, 'stable');

        % define the minimum number of trials each neuron must have for
        % each condition
        n_trials2use = 20;

        n_trials_per_cond_per_unit = data.n_trials_per_choice_cond;

        % find neurons with less than n_trials2use for any condition
        units2reject = find(sum(n_trials_per_cond_per_unit < n_trials2use) > 0);

        % identify the trials associated with these to-be-rejected units
        reject_ix  = zeros(numel(data.u_names), 1);
        for ix = 1:numel(units2reject)

            reject_ix(contains(data.u_names, neuron_names(units2reject(ix)))) = 1;

        end % of looping over neurons to reject

        % convert reject_ix to a boolean
        reject_ix = logical(reject_ix);

        % remove the trials associated with the rejected units
        bhv(reject_ix, :) = [];
        mean_fr(reject_ix) = [];
        unit_per_trial(reject_ix) = [];
        neuron_names(units2reject) = [];
        n_trials_per_cond_per_unit(:, units2reject) = [];

        boot_acc=NaN(2,2, n_boots);


        % initialize an array that will be used for this round's
        % decoder
        fr_matrix = NaN(n_trials2use*length(choice_trial_types), length(neuron_names), n_boots);

        % loop over each neuron
        pw = PoolWaitbar(numel(neuron_names), 'finding trials for each neuron...');
        for u = 1:numel(neuron_names)
        increment(pw);

            u_ix = contains(unit_per_trial, neuron_names(u));

            trials2keep = [];
            decoder_labels = [];

            % now loop over each of the choice_trial_types and find
            % n_trials2use random trials of this type for this neuron
            for c = 1:length(choice_trial_types)

                state_ix = bhv.state == choice_trial_types(c, 1);
                val_ix = bhv.chosenval == choice_trial_types(c, 2);

                c_trials = find(state_ix & val_ix & u_ix);
                n_c_trials = numel(c_trials);

                % shuffle c_trials 1000 times (for the bootstraps)
                shuff_c_trials = c_trials(randi(n_c_trials, [n_c_trials, n_boots]));

                trials2keep = [trials2keep ; shuff_c_trials(1:n_trials2use, :)];

                % create the labels for this condition
                state_label = ones(n_trials2use, 1).*choice_trial_types(c, 1);
                val_label = ones(n_trials2use, 1).*choice_trial_types(c, 2);
                state_val_label = state_label.*val_label;

                decoder_labels = [decoder_labels; state_label, val_label, state_val_label];

            end % of looping over conditions

            % accumulate the firing rates for each neuron in each condition
            % and bootstrap
            fr_matrix(:, u, :) = mean_fr(trials2keep);

        end % of looping over each neuron
        delete(pw);

%-----------------------------------------------------------------------------        

        % pre-determine the rows to use for the train/test split on each bootstrap 

        values_to_find = unique(decoder_labels(:,3), 'stable');
        % Find indices of the first occurrence of each value
        for i = 1:length(values_to_find)
            first_trial_per_class(i) = find(decoder_labels(:,3) == values_to_find(i), 1);
        end

        % now loop over number of neurons to construct a plot where the
        % x-axis is n_units and the y-axis is decoder_accuracy
        n_neurons2use = 5:5:numel(neuron_names);

        % pre-allocate arrays that get populated in the subsequent
        % parallelized bootstrapping
        cross_acc = NaN(2, 2, n_boots, numel(n_neurons2use));
        within_acc = NaN(n_boots, numel(n_neurons2use));
        between_acc = NaN(n_boots, numel(n_neurons2use));

        confusion_acc = NaN(8, 8, n_boots, numel(n_neurons2use));

        pw = PoolWaitbar(numel(n_neurons2use), 'Looping over n units to use...');
        parfor n_units2use = 1:numel(n_neurons2use)
        increment(pw);
            % generate n_boot iterations of which neurons to use
            units2use = randi(numel(neuron_names), [n_neurons2use(n_units2use), n_boots]);

            [cross_acc(:,:,:, n_units2use), within_acc(:,n_units2use), between_acc(:,n_units2use), confusion_acc(:, :, :, n_units2use)] =...
                cross_condition_generalization_decoder(units2use, n_boots, fr_matrix, n_trials2use, decoder_labels);
        end
        delete(pw);

        % store the results
        gcca_results.(subj_names{s}).(brain_areas{a}).all_acc = cross_acc;
        gcca_results.(subj_names{s}).(brain_areas{a}).within_acc = within_acc;
        gcca_results.(subj_names{s}).(brain_areas{a}).between_acc = between_acc;
        gcca_results.(subj_names{s}).(brain_areas{a}).confusion_acc = confusion_acc;

    end % of looping over brain areas

end % of looping over subjects

% save the chance level
gcca_results.chance_level = 1/numel(unique(choice_trial_types(:,2)));


% do some plotting
[D_OFC_within_mean, D_OFC_within_ci] = GetMeanCI(gcca_results.Don.OFC.within_acc, 'bootstrap');
[D_OFC_between_mean, D_OFC_between_ci] = GetMeanCI(gcca_results.Don.OFC.between_acc, 'bootstrap');

[D_HPC_within_mean, D_HPC_within_ci] = GetMeanCI(gcca_results.Don.HPC.within_acc, 'bootstrap');
[D_HPC_between_mean, D_HPC_between_ci] = GetMeanCI(gcca_results.Don.HPC.between_acc, 'bootstrap');


[K_OFC_within_mean, K_OFC_within_ci] = GetMeanCI(gcca_results.King.OFC.within_acc, 'bootstrap');
[K_OFC_between_mean, K_OFC_between_ci] = GetMeanCI(gcca_results.King.OFC.between_acc, 'bootstrap');

[K_HPC_within_mean, K_HPC_within_ci] = GetMeanCI(gcca_results.King.HPC.within_acc, 'bootstrap');
[K_HPC_between_mean, K_HPC_between_ci] = GetMeanCI(gcca_results.King.HPC.between_acc, 'bootstrap');


cmap = cbrewer('qual', 'Set1', 9);

figure; 
subplot(2,2,1)
hold on
shadedErrorBar([1:numel(D_OFC_within_mean)]*5, D_OFC_within_mean, D_OFC_within_ci, 'lineprops',{'color', cmap(2, :), 'LineWidth', 2});
shadedErrorBar([1:numel(D_HPC_within_mean)]*5, D_HPC_between_mean, D_HPC_within_ci, 'lineprops',{'color', cmap(9, :), 'LineWidth', 2});
xlabel('n units')
ylabel({'Subject D', 'Decoder Accuracy'})
ylim([.25, .6])
axis square

subplot(2,2,2)
hold on
shadedErrorBar([1:numel(D_OFC_between_mean)]*5, D_OFC_between_mean, D_OFC_between_ci, 'lineprops',{'color', cmap(2, :), 'LineWidth', 2});
shadedErrorBar([1:numel(D_HPC_between_mean)]*5, D_HPC_between_mean, D_HPC_between_ci, 'lineprops',{'color', cmap(9, :), 'LineWidth', 2});
xlabel('n units')
ylim([.25, .32])
axis square



subplot(2,2,3)
hold on
shadedErrorBar([1:numel(K_OFC_within_mean)]*5, K_OFC_within_mean, K_OFC_within_ci, 'lineprops',{'color', cmap(2, :), 'LineWidth', 2});
shadedErrorBar([1:numel(K_HPC_within_mean)]*5, K_HPC_between_mean, K_HPC_within_ci, 'lineprops',{'color', cmap(9, :), 'LineWidth', 2});
ylim([.25, .55])
ylabel({'Subject K', 'Decoder Accuracy'})
xlabel('n units')
axis square


subplot(2,2,4)
hold on
shadedErrorBar([1:numel(K_OFC_between_mean)]*5, K_OFC_between_mean, K_OFC_between_ci, 'lineprops',{'color', cmap(2, :), 'LineWidth', 2});
shadedErrorBar([1:numel(K_HPC_between_mean)]*5, K_HPC_between_mean, K_HPC_between_ci, 'lineprops',{'color', cmap(9, :), 'LineWidth', 2});
ylim([.25, .32])
xlabel('n units')
axis square




