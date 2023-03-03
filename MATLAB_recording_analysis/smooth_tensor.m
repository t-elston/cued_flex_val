function [out_data] = smooth_tensor(indata, win_size, dim)

[n_trials, n_times, n_units] = size(indata);

    choice_raw_firing_rates = smoothdata(choice_spk_table,'gaussian','omitnan','window',100) * 1000;




end % of function 