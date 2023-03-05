function [cue_amp, cue_phase, choice_amp, choice_phase, lfp_ts] = chop_lfp_v01(lfp,fs,freq_range, align_times, trials2use, offsets, ch_ix)

%-------------------------------------------------------------------------
% wrapper function to run chop_lfp_v01 in parallel
% Thomas Elston
% 14 Sept 2022
%-------------------------------------------------------------------------




 
    freqs = 2:50;
    n_freqs = numel(freqs);
    
    for ch_ix = 1:numel(channels_to_check)
        ch = channels_to_check(ch_ix);
        
        [fs, n, ~, ~, lfp] = plx_ad(OpenedFileName, names{ch});
        
        for f = 1:numel(freqs)
            
            disp([ 'ch: ' num2str(ch_ix) ' / ' num2str(numel(channels_to_check)) ', freq(Hz): ' num2str(freqs(f))])
            
            [cue_amp(f,:,:,ch_ix), cue_phase(f,:,:,ch_ix), choice_lfp(f,:,:,ch_ix), choice_phase(f,:,:,ch_ix), lfp_ts(f,:,:,ch_ix)] = ...
                chop_lfp_v01(lfp,fs,[freqs(f)-1,freqs(f)+1],[trialinfo(:,5),trialinfo(:,4)], trialinfo(:,3), [psth_neg_offset,psth_pos_offset]);
        end
        
    end % of looping over channels and extracting LF



end % of function