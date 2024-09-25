% Define parameters

ns_ofc_cue_A_z = cue_zval(ofc_ix & ~cue_sel ,1); 
ns_ofc_cue_B_z = cue_zval(ofc_ix & ~cue_sel ,2); 
ns_ofc_cue_A_z(~ns_choice_sel) = [];
ns_ofc_cue_B_z(~ns_choice_sel) = [];

n = numel(ns_ofc_cue_A_z); % Number of items

% Define the threshold for determining selectivity
threshold = 2;

% Calculate the occurrences of selective neurons in each condition
selective_in_A = (ns_ofc_cue_A_z > threshold) & (ns_ofc_cue_B_z < threshold);
selective_in_B = (ns_ofc_cue_A_z < threshold) & (ns_ofc_cue_B_z > threshold);
selective_in_both = (ns_ofc_cue_A_z > threshold) & (ns_ofc_cue_B_z > threshold);

% Calculate the proportion of neurons selective in both conditions
proportion_selective_in_both = sum(selective_in_both) / n;

% Calculate the union of the two probabilities
expected_proportion_selective_in_both = mean(selective_in_A) + mean(selective_in_B) - proportion_selective_in_both;

% Calculate the actual proportion of neurons selective in both conditions
actual_proportion_selective_in_both = sum(selective_in_both) / n;

% Perform binomial test to compare the observed proportion to the expected chance level
p_value = 1 - binocdf(16, 113, expected_proportion_selective_in_both, 'upper');

% Display the results
fprintf('Expected proportion of neurons selective in both conditions by chance: %.4f\n', expected_proportion_selective_in_both);
fprintf('Actual proportion of neurons selective in both conditions: %.4f\n', actual_proportion_selective_in_both);
fprintf('Binomial test p-value: %.4f\n', p_value);

% Calculate the occurrences of neurons selective in either condition
selective_in_either = sum((z_stat_condition1 > threshold) | (z_stat_condition2 > threshold));