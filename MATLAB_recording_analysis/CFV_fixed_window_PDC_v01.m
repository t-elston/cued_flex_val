function [outPDC] = CFV_fixed_window_PDC_v01(OFC_trial_lfp, HPC_trial_lfp, FOI, fs)


% convert FOI to sample intervals
freqs = round(fs./FOI);

[n_times, n_trials] = size(OFC_trial_lfp);

hpc2ofc = NaN(n_trials, numel(freqs));
ofc2hpc = NaN(n_trials, numel(freqs));


for trial = 1:n_trials
      
    % cat two waves into one array
    waveboth= cat(2,OFC_trial_lfp(:, trial),HPC_trial_lfp(:, trial));   
    
    shuffle_waveboth = cat(2,shuffle(OFC_trial_lfp(:, trial)),shuffle(HPC_trial_lfp(:, trial)));   
    
    [pdc] = PDC(waveboth,15,freqs);
    shuffle_pdc = PDC(shuffle_waveboth,15,freqs);
    
    hpc2ofc(trial, :)= pdc(:,1,2);
    ofc2hpc(trial, :)= pdc(:,2,1);
    
    s_hpc2ofc(trial, :)= shuffle_pdc(:,1,2);
    s_ofc2hpc(trial, :)= shuffle_pdc(:,2,1);

    
end % of looping over trials

outPDC(:,1) = smoothdata(circshift(nanmean(ofc2hpc),2),'gaussian',1);
outPDC(:,2) = smoothdata(circshift(nanmean(hpc2ofc),2),'gaussian',1);
outPDC(:,3) = smoothdata(circshift(nanmean(s_ofc2hpc),2),'gaussian',1);
outPDC(:,4) = smoothdata(circshift(nanmean(s_hpc2ofc),2),'gaussian',1);
outPDC(:,5) = outPDC(:,2) - outPDC(:,1);
outPDC(:,6) = outPDC(:,4) - outPDC(:,3);


end % of function 