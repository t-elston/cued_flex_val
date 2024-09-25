function [cross_acc, within_acc, between_acc, confusion_acc] = cross_condition_generalization_decoder(units2use, n_boots, fr_matrix, n_trials2use, decoder_labels)


cross_acc = NaN(2,2, n_boots);
within_acc = NaN(n_boots, 1);
between_acc = NaN(n_boots, 1);

confusion_acc = NaN(8, 8, n_boots);

% train on state 1, test on state -1
pos1_ix = decoder_labels(:,1) == 1;
neg1_ix = decoder_labels(:,1) == -1;

% pre-determine the rows to use for the train/test split on each bootstrap

values_to_find = unique(decoder_labels(:,3), 'stable');
% Find indices of the first occurrence of each value
for i = 1:length(values_to_find)
    first_trial_per_class(i) = find(decoder_labels(:,3) == values_to_find(i), 1);
end


% loop over the bootstraps for this number of neurons
for b = 1:n_boots

    b_fr = fr_matrix(:, units2use(:, b), b);

    [~, pcs, ~, ~, explained_var] = pca(b_fr);

    % find the PCs that explain 95% of the variance
    cum_exp_var = cumsum(explained_var);
    pc95 = min(find(cum_exp_var >= 95));

    % keep only those most informative PCs
    pcs = pcs(:, 1: pc95);

    pcs = b_fr;
    

    % use a leave one out train/test testing procedure
    acc = NaN(2,2,n_trials2use);
    b_confusion_acc = NaN(8,8,n_trials2use); 

    for f = 1:n_trials2use
        % define the testing set
        test_ix = zeros(numel(pos1_ix), 1);
        test_ix(first_trial_per_class + f - 1) = 1;

        % TRAIN/TEST within/across individual states
        % train the SVM decoders
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

        % TRAIN/TEST an 8-way decoder
        full_decoder = fitcecoc(pcs(~test_ix, :), decoder_labels(~test_ix, 3));
        predicted_labels = predict(full_decoder, pcs(logical(test_ix), :));
        [C, order] = confusionmat(decoder_labels(logical(test_ix), 3), predicted_labels,...
            'Order',[1 2 3 4, -1, -2, -3, -4]);

        b_confusion_acc(:,:, f) = C;


    end % loop over the folds of a leave-one-out-approach

    b_acc = mean(acc,3);
    confusion_acc(:,:, b) = mean(b_confusion_acc, 3);

    cross_acc(:,:, b) = b_acc;
    within_acc(b, 1) = mean([b_acc(1,1), b_acc(2,2)]);
    between_acc(b, 1) = mean([b_acc(1,2), b_acc(2,1)]);


end % of looping over bootstraps

xx=[];



