% Choice period decoding analysis

base_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\preprocessed data/';

subj_names = {'Don', 'King'};
brain_areas = {'HPC', 'OFC'};


% details of an 8 way decoder
choice_trial_types = [ 1, 1, 1, 1,-1,-1,-1,-1; % states
                       1, 2, 3, 4, 1, 2, 3, 4]'; % values

% how many times should we run the decoder?
n_boots = 10;

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
            [~, win_end_ix] = min(abs(data.ts - 500));

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


            % prepare to loop over different numbers of randomly
            % selected neurons
            n_units = numel(neuron_names);

            neuron_increments = 5:5:n_units;

            within_acc=NaN(n_boots, numel(neuron_increments));
            between_acc=NaN(n_boots, numel(neuron_increments));

            % start looping over bootstraps
            pw = PoolWaitbar(n_boots, 'bootstrapping...');
            for b = 1:n_boots
                increment(pw);

                % initialize an array that will be used for this round's
                % decoder
                fr_matrix = NaN(n_trials2use*length(choice_trial_types), length(neuron_names));

                % loop over each neuron
                for u = 1:numel(neuron_names)

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

                        % shuffle c_trials and grab a random n_trials2use
                        shuff_c_trials = c_trials(randi(n_c_trials, [n_c_trials, 1]));

                        trials2keep = [trials2keep ; shuff_c_trials(1:n_trials2use)];

                        % create the labels for this condition
                        state_label = ones(n_trials2use, 1).*choice_trial_types(c, 1);
                        val_label = ones(n_trials2use, 1).*choice_trial_types(c, 2);
                        state_val_label = state_label.*val_label;

                        decoder_labels = [decoder_labels; state_label, val_label, state_val_label];

                    end % of looping over conditions

                    % add these selected trials to the big array of firing rates
                    fr_matrix(:, u) = mean_fr(trials2keep);

                end % of looping over each neuron
                
                % zscore the firing rates for this run
                fr_matrix = zscore(fr_matrix, [], 1);


                % now that we've accumulated the firing rates for all
                % neurons on this bootstrap, let's run the decoder
                % train on state 1, test on state -1
                pos1_ix = decoder_labels(:,1) == 1;
                neg1_ix = decoder_labels(:,1) == -1;

                values_to_find = unique(decoder_labels(:,3), 'stable');
                % Find indices of the first occurrence of each value
                for i = 1:length(values_to_find)
                    first_trial_per_class(i) = find(decoder_labels(:,3) == values_to_find(i), 1);
                end

                % now loop over neurons
                for u_ix = 1:numel(neuron_increments)

                    units2use = randi(n_units, [neuron_increments(u_ix), 1]);

                    % select the units for this round
                    b_fr = fr_matrix(:, units2use);

                    % reduce the dimensionality with PCA
                    [coeff, pcs, latent, ~, explained_var] = pca(b_fr);

                    % find the PCs that explain 90% of the variance
                    cum_exp_var = cumsum(explained_var);
                    pc90 = min(find(cum_exp_var >= 90));

                    % keep only those most informative PCs
                    pcs = pcs(:, 1: pc90);

                    % use a leave one out train/test testing procedure
                    acc = NaN(2,2,n_trials2use);
                    for f = 1:n_trials2use
                        % define the testing set
                        test_ix = zeros(numel(pos1_ix), 1);
                        test_ix(first_trial_per_class + f - 1) = 1;

                        % train the decoders
                        % pos1_decoder = fitcdiscr(pcs(pos1_ix & ~test_ix, :),...
                        %     decoder_labels(pos1_ix & ~test_ix, 2));
                        % neg1_decoder = fitcdiscr(pcs(neg1_ix & ~test_ix, :),...
                        %     decoder_labels(neg1_ix & ~test_ix, 2));

                        pos1_decoder = fitcecoc(pcs(pos1_ix & ~test_ix, :),...
                            decoder_labels(pos1_ix & ~test_ix, 2));
                        neg1_decoder = fitcecoc(pcs(neg1_ix & ~test_ix, :),...
                            decoder_labels(neg1_ix & ~test_ix, 2));

                        % test the decoders
                        acc(1,1,f) = mean(predict(pos1_decoder,...
                            pcs(pos1_ix & test_ix, :)) == decoder_labels(pos1_ix & test_ix, 2));

                        acc(1,2,f) = mean(predict(pos1_decoder,...
                            pcs(neg1_ix & test_ix, :)) == decoder_labels(neg1_ix & test_ix, 2));

                        acc(2,1,f) = mean(predict(neg1_decoder,...
                            pcs(pos1_ix & test_ix, :)) == decoder_labels(pos1_ix & test_ix, 2));

                        acc(2,2,f) = mean(predict(neg1_decoder,...
                            pcs(neg1_ix & test_ix, :)) == decoder_labels(neg1_ix & test_ix, 2));
                    end % loop over the folds of a leave-one-out-approach

                    boot_acc = mean(acc,3);

                    within_acc(b, u_ix) = nanmean([boot_acc(1,1), boot_acc(2,2)]);
                    between_acc(b, u_ix) = nanmean([boot_acc(1,2), boot_acc(2,1)]);

                end % of looping over this many neurons for this bootstrap

            end % of looping over bootstraps
            delete(pw);

            % store the results
            gcca_results.(subj_names{s}).(brain_areas{a}).within_acc = within_acc;
            gcca_results.(subj_names{s}).(brain_areas{a}).across_acc = between_acc;

        end % of looping over brain areas

end % of looping over subjects

% save the chance level
gcca_results.chance_level = 1/numel(unique(choice_trial_types(:,2)));
