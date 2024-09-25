function [GC_out] = freq_resolved_GC(inwave1, inwave2, fs, align_events)

freq_bands = [ 2, 8;
              10, 15;
              15, 30;
              40, 70];

% loop over the frequencies
for f = 1:size(freq_bands, 1)

    f_hpc = bandpass_filter(inwave1, fs, freq_bands(f, 1), freq_bands(f, 2), 4);
    f_ofc = bandpass_filter(inwave2, fs, freq_bands(f, 1), freq_bands(f, 2), 4);

    [OFC_trial_lfp, lfp_ts] = chop_raw_lfp_v01(f_ofc, fs, align_events, [0, 1000]);
    [HPC_trial_lfp, ~] = chop_raw_lfp_v01(f_hpc, fs, align_events, [0, 1000]);

    OFC_trial_lfp = zscore(OFC_trial_lfp);
    HPC_trial_lfp = zscore(HPC_trial_lfp);

    % now compute the granger causality for each trial
    for t = 1:size(HPC_trial_lfp, 2)
        
        [GC_xy(t, f), GC_yx(t, f)] = GrangerCausality(HPC_trial_lfp(:, t), OFC_trial_lfp(:, t), 15);

        % now shuffle the waves
        [s_GC_xy(t, f), s_GC_yx(t, f)] = GrangerCausality(shuffle(HPC_trial_lfp(:, t)), shuffle(OFC_trial_lfp(:, t)), 15);

    end % of looping over trials

end % of looping over frequencies

GC_out(:, 1) = mean(GC_yx - GC_xy)';
GC_out(:, 2) = mean(s_GC_yx - s_GC_xy)';


end % of function