%%% aaa_psuedopop_analysis_v01

% first, aggregate the data across the files for the pseudo-pop
datadir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\recording_data\';

n_boots = 100;
[train_labels, train_firing_rates, test_labels, test_firing_rates, OFC_ix] = pseudo_aggregate_data_v01(n_boots, datadir);

[HPC_accs, OFC_accs] = pseudo_LDA_v01(train_labels, train_firing_rates, test_labels, test_firing_rates, OFC_ix, n_boots);