function ch_details = compute_channel_PSDs_v3(ch_details, pl2_fname)

n_ch = numel(ch_details.Tower);

    % FOOOF settings
    settings = struct();  % Use defaults
    f_range = [1, 50];

for ch_ix = 1:n_ch
    
    % compute PSD
    ChLFP = PL2Ad(pl2_fname, ch_details.ch_name{ch_ix});
    
    [psd, freqs] = pwelch(ChLFP.Values, 1000, [], [1:50], 1000);

    % Fit a power-law model (1/f) using linear regression in log-log space
    log_freqs = log10(freqs);
    log_psd = log10(psd);
    coefficients = polyfit(log_freqs, log_psd, 1);
    fit_line = polyval(coefficients, log_freqs);

    % Subtract the fitted line from the original PSD to remove the aperiodic component
    psd_detrended = psd - 10.^(fit_line);

    ch_details.theta_power(ch_ix) = mean(psd_detrended(4:8));
    ch_details.psd(ch_ix, :) = psd_detrended;
     
    ChLFP=[];
end % of looping over channels and finding which have theta



end % of function