function [boot_mean, boot_CI] = get_shuffled_confidence_intervals_for_coherence(indata, n_boots)

% shuffle data along the temporal dimension (dim = 2)

% get size of indata
[n_pairs, n_times] = size(indata); 

% create random sets of indices for the shuffle
shuff_indices = randi(n_times, n_boots, n_times);

b_mean = NaN(n_boots, n_times);

% loop over the bootstraps
for b = 1:n_boots
    
    b_mean(b, :) = mean(indata(:, shuff_indices(b,:)));

end % of looping over bootstraps

boot_mean = mean(b_mean);
boot_CI = std(b_mean);

end % of function