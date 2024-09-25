% Choice period decoding analysis

base_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\preprocessed data/';

subj_names = {'Don', 'King'};
brain_areas = {'HPC', 'OFC'};


% details of a 2 way state decoder
choice_trial_types = [ 1, 1, -1, -1; % states
                       1, -1, 1, -1]'; % types

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

        % get the mean firing rate over the cue period
        [~, win_start_ix] = min(abs(data.ts - -750));
        [~, win_end_ix] = min(abs(data.ts - -200));

        % mean firing rate for every trial of every neuron
        mean_fr = nanmean(data.zFRs(:, win_start_ix: win_end_ix), 2);

        % get the names of the neurons
        unit_per_trial = data.u_names;
        neuron_names = unique(unit_per_trial, 'stable');

        % define the minimum number of trials each neuron must have for
        % each condition
        n_trials2use = 200;

        boot_acc=NaN(n_boots, 2);

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
                type_ix = bhv.state_type == choice_trial_types(c, 2);

                c_trials = find(state_ix & type_ix & u_ix);
                n_c_trials = numel(c_trials);

                % shuffle c_trials 1000 times (for the bootstraps)
                shuff_c_trials = c_trials(randi(n_c_trials, [n_c_trials, n_boots]));

                trials2keep = [trials2keep ; shuff_c_trials(1:n_trials2use, :)];

                % create the labels for this condition
                state_label = ones(n_trials2use, 1).*choice_trial_types(c, 1);
                type_label = ones(n_trials2use, 1).*choice_trial_types(c, 2);
                state_type_labels = ones(n_trials2use, 1).*c;

                decoder_labels = [decoder_labels; state_label, type_label, state_type_labels];

            end % of looping over conditions

            % accumulate the firing rates for each neuron in each condition
            % and bootstrap
            fr_matrix(:, u, :) = mean_fr(trials2keep);

        end % of looping over each neuron
        delete(pw);

%-----------------------------------------------------------------------------        

        % now loop over number of neurons to construct a plot where the
        % x-axis is n_units and the y-axis is decoder_accuracy
        n_neurons2use = 50:50:numel(neuron_names);
        n_neurons2use = numel(neuron_names);

        % pre-allocate arrays that get populated in the subsequent
        % parallelized bootstrapping
        state_acc = NaN(n_boots, numel(n_neurons2use));

        pw = PoolWaitbar(numel(n_neurons2use), 'Looping over n units to use...');
        parfor n_units2use = 1:numel(n_neurons2use)
        increment(pw);

            % generate n_boot iterations of which neurons to use
            units2use = randi(numel(neuron_names), [n_neurons2use(n_units2use), n_boots]);

            state_acc(:, n_units2use) = parallel_state_decoder(units2use, n_boots, fr_matrix, n_trials2use, decoder_labels);
        
        end
        delete(pw);

        % store the results
        gcca_results.(subj_names{s}).(brain_areas{a}).state_acc = state_acc;
        gcca_results.(subj_names{s}).(brain_areas{a}).n_units = n_neurons2use;

    end % of looping over brain areas

end % of looping over subjects

cmap = cbrewer('qual', 'Set1', 9);




