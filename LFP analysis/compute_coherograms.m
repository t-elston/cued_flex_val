function [trial_coh, z_coh, coh_ts, outPDC, GC_out] = compute_coherograms(ch_tbl, pl2_fname, pl2_event_times)

%-----------------------------------------
% PARAMETERS for multi taper coherogram
fs=1000; % sampling freq
window=fs;
tapers=4;
detrend=1;
keep=3;
nFFT=window;
% noverlap = round(.9*window);
noverlap = round(window - (100/window)*window);
FOI = [1: 30]; % freq range of interest
%-----------------------------------------

% make filter for bandpassing LFP into range of 1 - 100 Hz
lfp_filt = designfilt('bandpassfir', ... response type
    'FilterOrder',1000, ... filter order
    'CutoffFrequency1', 1,... lower range
    'CutoffFrequency2', 100,... upper range
    'SampleRate', 1000); ... sampling rate


OFC_ch = ch_tbl(contains(ch_tbl.Target,'OFC'),:);
HPC_ch = ch_tbl(contains(ch_tbl.Target,'HPC'),:);

% make a table of the channel pairs
pair_tbl = table;
pair_ctr=0;
for ofc_ch_ix = 1:numel(OFC_ch.Tower)  
    for hpc_ch_ix = 1:numel(HPC_ch.Tower)
        
        pair_ctr = pair_ctr+1;
        
        pair_tbl.OFC_ch(pair_ctr) = OFC_ch.ch_name(ofc_ch_ix);
        pair_tbl.HPC_ch(pair_ctr) = HPC_ch.ch_name(hpc_ch_ix);
        
    end % of looping over HPC channels    
end % of looping over OFC channels


trial_coh={};
coh_ts = {};
z_coh = {};
outPDC=[];
GC_out=[];

if ~isempty(pair_tbl)
    
    n_pairs = numel(pair_tbl.OFC_ch);

    % now loop over the channel pairs
    trial_coh=cell(n_pairs,1);
    coh_ts= cell(n_pairs,1);
    outPDC = NaN(30, 6, n_pairs);
    GC_out = NaN(4, 2, n_pairs);
    pw = PoolWaitbar(n_pairs, 'Assessing OFC-HPC coherence...');
    parfor (pair_ix = 1:n_pairs, 2)    
    % for pair_ix = 1:n_pairs    

        increment(pw);
         
        OFC_lfp = PL2Ad(pl2_fname, pair_tbl.OFC_ch{pair_ix});
        HPC_lfp = PL2Ad(pl2_fname, pair_tbl.HPC_ch{pair_ix});
        fs = OFC_lfp.ADFreq;
          
        OFC_signal = fftfilt(lfp_filt.Coefficients,OFC_lfp.Values);
        HPC_signal = fftfilt(lfp_filt.Coefficients,HPC_lfp.Values);
        
        lag = 1+(1+1)*OFC_lfp.ADFreq/2;
        
        OFC_signal = OFC_signal(lag:end);
        HPC_signal = HPC_signal(lag:end);
         
        % chop the LFP into trials for PDC analysis
        [OFC_trial_lfp, lfp_ts] = chop_raw_lfp_v01(OFC_signal, OFC_lfp.ADFreq, pl2_event_times.cue_on, [0, 1000]); % was 1400
        [HPC_trial_lfp, ~] = chop_raw_lfp_v01(HPC_signal, HPC_lfp.ADFreq, pl2_event_times.cue_on, [0, 1000]);
        
        % now do a PDC analysis
        [outPDC(:,:,pair_ix)] = CFV_fixed_window_PDC_v01(OFC_trial_lfp, HPC_trial_lfp, FOI, OFC_lfp.ADFreq);  

        % now do freq-resolved Granger Causality
        [GC_out(:,:, pair_ix)] = freq_resolved_GC(HPC_signal, OFC_signal, fs, pl2_event_times.cue_on);
        
        % compute the cohereogram
        [yo, freqs, ts]=mtchg_vTE({[OFC_lfp.Values, HPC_lfp.Values],nFFT,fs,window,noverlap,tapers,detrend,keep},FOI);
    
        % extract spectra of interest from 4-D array
        mtCoh=yo(:,:,1,2);
    
        % chop coherogram into trials
        [trial_coh{pair_ix,1}, coh_ts{pair_ix,1}, z_coh{pair_ix,1}] = chop_coherence_into_trials_v02(mtCoh,round(ts*1000), pl2_event_times, freqs);
    
    end % of looping over channel pairs
    delete(pw);
    coh_ts = coh_ts{1};

end % of checking if the channel pair table was empty



        
end % of function
