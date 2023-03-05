function [trial_data, iti_data, coh_ts, outPDC, mean_lag] = CFV_compute_coherogram(ch_tbl, file_ch_names, full_pl2_fname, session_name, trialinfo)

%-----------------------------------------
% PARAMETERS for multi taper coherogram
fs=1000; % sampling freq
window=1*fs;
tapers=4;
detrend=1;
keep=3;
nFFT=window;
noverlap = round(.9*fs);
FOI = [1: 50]; % freq range of interest
%-----------------------------------------

% make filter for bandpassing LFP into theta range
lfp_filt = designfilt('bandpassfir', ... response type
    'FilterOrder',1000, ... filter order
    'CutoffFrequency1', 1,... lower range
    'CutoffFrequency2', 100,... upper range
    'SampleRate', 1000); ... sampling rate



session_channels = ch_tbl(contains(ch_tbl.file,session_name),:);
OFC_ch = session_channels(contains(session_channels.brain_area,'OFC'),:);
HPC_ch = session_channels(contains(session_channels.brain_area,'HPC'),:);

% make a table of the channel pairs
pair_tbl = table;
pair_ctr=0;
for ofc_ch_ix = 1:numel(OFC_ch.ch_num)
    
    for hpc_ch_ix = 1:numel(HPC_ch.ch_num)
        
        pair_ctr = pair_ctr+1;
        
        pair_tbl.OFC_ch(pair_ctr) = OFC_ch.ch_name(ofc_ch_ix);
        pair_tbl.HPC_ch(pair_ctr) = HPC_ch.ch_name(hpc_ch_ix);
        
    end % of looping over HPC channels
    
end % of looping over OFC channels


trial_data={};
iti_data = {};
coh_ts = {};
outPDC=[];

if ~isempty(pair_tbl)
    
    n_pairs = numel(pair_tbl.OFC_ch);

    % now loop over the channel pairs
    trial_data=cell(n_pairs,1);
    iti_data = cell(n_pairs,1);
    coh_ts= cell(n_pairs,1);
    outPDC = NaN(50,6, n_pairs);
    pw = PoolWaitbar(n_pairs, 'Assessing OFC-HPC coherence...');
    parfor (pair_ix = 1:n_pairs, 2)    
%     for pair_ix = 1:n_pairs    

         increment(pw);
         
        OFC_lfp = PL2Ad(full_pl2_fname, pair_tbl.OFC_ch{pair_ix});
        HPC_lfp = PL2Ad(full_pl2_fname, pair_tbl.HPC_ch{pair_ix});
        
        OFC_signal = OFC_lfp.Values;
        HPC_signal = HPC_lfp.Values;
        
        OFC_signal = fftfilt(lfp_filt.Coefficients,OFC_lfp.Values);
        HPC_signal = fftfilt(lfp_filt.Coefficients,HPC_lfp.Values);
        
        lag = 1+(1+1)*OFC_lfp.ADFreq/2;
        
        OFC_signal = OFC_signal(lag:end);
        HPC_signal = HPC_signal(lag:end);
        


        
        % chop the LFP into trials for PDC analysis
        [OFC_trial_lfp, lfp_ts] = chop_raw_lfp_v01(OFC_signal, OFC_lfp.ADFreq, trialinfo(:,5), [0, 1400]);
        [HPC_trial_lfp, ~] = chop_raw_lfp_v01(HPC_signal, HPC_lfp.ADFreq, trialinfo(:,5), [0, 1400]);

        % now do a sliding-window PDC analysis
        [outPDC(:,:,pair_ix)] = CFV_fixed_window_PDC_v01(OFC_trial_lfp, HPC_trial_lfp, [1:50], OFC_lfp.ADFreq);        
    
%         % compute the cohereogram
        [yo,freqs,ts]=mtchg_vTE({[OFC_lfp.Values, HPC_lfp.Values] ,nFFT,fs,window,noverlap,tapers,detrend,keep},FOI);
    
        %extract spectra of interest from 4-D array
        mtCoh=yo(:,:,1,2);
    
        % chop coherogram into trials
        [trial_data{pair_ix,1}, iti_data{pair_ix,1}, coh_ts{pair_ix,1}] = chop_coherence_into_trials_v01(mtCoh,round(ts*1000),trialinfo);
    
    
    end % of looping over channel pairs
    delete(pw);
    coh_ts = coh_ts{1};

end



        
end % of function
