function [trial_coh, coh_ts, outPDC] = phase_based_coherence(channels2use, pl2_fname, pl2_event_times)

OFC_ch = channels2use(contains(channels2use.Target,'OFC'),:);
HPC_ch = channels2use(contains(channels2use.Target,'HPC'),:);

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


if ~isempty(pair_tbl)
    
    % create a bunch of filters to assess the data within
    filters = {};
    
    for f = 2:31
        
        l_bound = f-1;
        
        u_bound = f+1;
        
        filters{f} = designfilt('bandpassfir', ... response type
            'FilterOrder',1000, ... filter order
            'CutoffFrequency1', l_bound,... lower range
            'CutoffFrequency2', u_bound,... upper range
            'SampleRate', 1000); ... sampling rate
            
    end % of looping over filters to construct
    filters(1) = []; % remove the first element because it's empty
    
    
    
    
    trial_coh={};
    coh_ts = {};
    outPDC=[];
    
    

    n_pairs = numel(pair_tbl.OFC_ch);
    
    % now loop over the channel pairs
    trial_coh=cell(n_pairs,1);
    coh_ts= cell(n_pairs,1);
    outPDC = NaN(50,6, n_pairs);
    pw = PoolWaitbar(n_pairs, 'Assessing OFC-HPC coherence...');
    %     parfor (pair_ix = 1:n_pairs, 2)
    for pair_ix = 1:n_pairs
        
        increment(pw);
        
        OFC_lfp = PL2Ad(pl2_fname, pair_tbl.OFC_ch{pair_ix});
        HPC_lfp = PL2Ad(pl2_fname, pair_tbl.HPC_ch{pair_ix});
        
        for f = 1:numel(filters)
            
        % extract this filter
        fq_filt = filters{f};
        
        % apply the bandpass
        OFC_signal = fftfilt(fq_filt.Coefficients,OFC_lfp.Values);
        HPC_signal = fftfilt(fq_filt.Coefficients,HPC_lfp.Values);
        
        % adjust time to account for filter lag
        LFPtime = 1:numel(OFC_lfp.Values); %ms
        lag = 1+(1+1)*OFC_lfp.ADFreq/2;
        OFC_signal = OFC_signal(lag:end);
        HPC_signal = HPC_signal(lag:end);
        
        % extract the instantaneous phase
        ofc_phase = angle(hilbert(OFC_signal));
        hpc_phase = angle(hilbert(HPC_signal));
        
        % chop the LFP phase into trials
        [ofc_trials] = chop_lfp_into_trials(ofc_phase, pl2_event_times, 2000);
        [hpc_trials] = chop_lfp_into_trials(hpc_phase, pl2_event_times, 2000);
        
        % assess the phase-based coherence
        phase_coh(f, :, pair_ix) = compute_LFP_phase_alignment(ofc_trials - hpc_trials);
        
%         % now do a PDC analysis
%         [outPDC(:,:,pair_ix)] = CFV_fixed_window_PDC_v01(OFC_trial_lfp, HPC_trial_lfp, [1:50], OFC_lfp.ADFreq);
        
        end % of looping over frequency bands
    end % of looping over channel pairs
    delete(pw);
    coh_ts = coh_ts{1};
    
end % of checking if the channel pair table was empty



        
end % of function
