function [outPDC] = CFV_sliding_window_PDC_v01(OFC_trial_lfp, HPC_trial_lfp, FOI, fs)

% convert FOI to sample intervals
freqs = round(fs./FOI);

[n_times, n_trials] = size(OFC_trial_lfp);

PDC_ts = [];

t_step = 100;

PDC_ts = 0:t_step:n_times;

hpc2ofc = NaN(numel(freqs), numel(PDC_ts), n_trials);
ofc2hpc = NaN(numel(freqs), numel(PDC_ts), n_trials);


for trial = 1:n_trials
    ctr=0;
    for time = 0:t_step:n_times
        ctr = ctr+1;
        % get the details of this trial window
        
        if time - 499 > 0      
            t_start = time-499;
        else
            t_start = 1;
        end
                    
        if time + 499 < n_times      
            t_end = time + 499;
        else
            t_end = n_times;
        end
        
        
        PDC_ts(ctr) = time;
        
        waveboth= cat(2,OFC_trial_lfp(t_start:t_end, trial),HPC_trial_lfp(t_start:t_end, trial));   % cat two waves into one array
        pdc=[];
        [pdc] = PDC(waveboth,20,freqs);
        OFC_to_HPC= pdc(:,1,2);       
        HPC_to_OFC= pdc(:,2,1);
        
        % now put these values into a trial array 
        hpc2ofc(:,ctr, trial) = HPC_to_OFC;
        ofc2hpc(:,ctr, trial) = OFC_to_HPC;


        
        
        
        
    end % of looping over times
    
end % of looping over times





end % of function 