% Define parameters

ns_ofc_cue_A_z = cue_zval(ofc_ix & ~cue_sel ,1); 
ns_ofc_cue_B_z = cue_zval(ofc_ix & ~cue_sel ,2); 
ns_ofc_cue_A_z(~ns_choice_sel) = [];
ns_ofc_cue_B_z(~ns_choice_sel) = [];

% Define Z-statistics for each condition
z_stat_condition1 = ns_ofc_cue_A_z; 
z_stat_condition2 = ns_ofc_cue_B_z; 
n = numel(z_stat_condition2); % Number of items

% Define the threshold for determining selectivity
threshold = 2;

% Calculate the occurrences of selective neurons in each condition
selective_in_1 = ((z_stat_condition1 > threshold) & (z_stat_condition2 < threshold)) | ((z_stat_condition1 < threshold) & (z_stat_condition2 > threshold));
selective_in_both = (z_stat_condition1 > threshold) & (z_stat_condition2 > threshold);

% Calculate the proportion of neurons selective in both conditions
proportion_selective_in_both = sum(selective_in_both) / n;

% Calculate the expected proportion of neurons selective in one condition but not in both conditions by chance
expected_proportion_selective_in_1 = 2 * (1 - proportion_selective_in_both) * proportion_selective_in_both;

% Calculate the actual proportion of neurons selective in one condition but not in both conditions
actual_proportion_selective_in_1 = sum(selective_in_1 & ~selective_in_both) / n;

% Perform binomial test to compare the observed proportion to the expected chance level
p_value = binocdf(sum(selective_in_1 & ~selective_in_both), n, expected_proportion_selective_in_1, 'upper');

% Display the results
fprintf('Expected proportion of neurons selective in one condition but not in both conditions by chance: %.4f\n', expected_proportion_selective_in_1);
fprintf('Actual proportion of neurons selective in one condition but not in both conditions: %.4f\n', actual_proportion_selective_in_1);
fprintf('Binomial test p-value: %.4f\n', p_value);

selective_in_either = sum((z_stat_condition1 > threshold) | (z_stat_condition2 > threshold));
