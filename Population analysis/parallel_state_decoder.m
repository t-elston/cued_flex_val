function [state_acc] = parallel_state_decoder(units2use, n_boots, fr_matrix, n_trials2use, decoder_labels)


state_acc = NaN(n_boots, 1);

% pre-determine the rows to use for the train/test split on each bootstrap
values_to_find = unique(decoder_labels(:,3), 'stable');
% Find indices of the first occurrence of each value
for i = 1:length(values_to_find)
    first_trial_per_class(i) = find(decoder_labels(:,3) == values_to_find(i), 1);
end

state_A_ix = decoder_labels(:, 1) == -1;
state_B_ix = decoder_labels(:, 1) == 1;



% loop over the bootstraps for this number of neurons
for b = 1:n_boots

    b_fr = fr_matrix(:, units2use(:, b), b);

    % [~, pcs, ~, ~, explained_var] = pca(b_fr);
    % % 
    % % % find the PCs that explain 95% of the variance
    % cum_exp_var = cumsum(explained_var);
    % pc95 = min(find(cum_exp_var >= 95));
    % % 
    % % % keep only those most informative PCs
    % pcs = pcs(:, 1: pc95);
    pcs = b_fr;
    

    % use a stratified leave one out train/test testing procedure
    acc = NaN(n_trials2use, 1);

    for f = 1:n_trials2use

        % define the testing set
        test_ix = zeros(size(pcs,1), 1);
        test_ix(first_trial_per_class + f - 1) = 1;
        test_ix = logical(test_ix);

        % TRAIN/TEST within/across individual states
        % train the SVM decoder
        svm = fitcdiscr(pcs(~test_ix, :), decoder_labels(~test_ix, 1));

        % test the decoder
        acc(f,1) = mean(predict(svm, pcs(test_ix, :)) == decoder_labels(test_ix, 1));


    end % loop over the folds of a leave-one-out-approach

    state_acc(b, 1) = mean(acc);

end % of looping over bootstraps

xx=[];



