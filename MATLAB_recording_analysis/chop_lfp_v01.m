function [cue_amp, cue_phase, choice_amp, choice_phase, lfp_ts] = chop_lfp_v01(lfp,fs,freq_range, align_times, trials2use, offsets)

%-------------------------------------------------------------------------
% chops lfp into trials based on the timestamps in align_times and
% trials2use
% Thomas Elston
% 14 Sept 2022
%-------------------------------------------------------------------------
% INPUTS
% lfp         - datastream of unfiltered local field potential
% fs          - sampling frequency of lfp (e.g. 1000 = 1000Hz)
% freq_range  - 2-element array where the first entry is the lower bound of
%               a bandpass and the 2nd entry is the upper bound (e.g. [2,4] with bandpass
%               between 2 and 4 Hz
% align_times - n_trials x n_align_events array of times to align data to organized such that each
%               row corresponds to a single trial and each column
%               corresponds to a different aligning event
% trials2use  - index of which trials to use in align_times
% offsets     - 2-element array detailing how many milliseconds 
%               [before,after] align_times to use when constructing the LFP scalogram

% OUTPUTS
% cue_amp      - an n_trials x n_times array aligned to cue
% choice_amp   - an n_trials x n_times array aligned to choice on
% lfp_ts       - an n_times array corresponding to the time relative to
%                the aligning event
%-------------------------------------------------------------------------

n_trials = numel(trials2use);

% compute time stamps
lfp_ts =  -1*offsets(1): offsets(2)-1;

% initialize output
cue_erp = NaN(numel(lfp_ts), n_trials);
choice_erp = NaN(numel(lfp_ts), n_trials);

cue_amp = NaN(numel(lfp_ts), n_trials);
choice_amp = NaN(numel(lfp_ts), n_trials);
cue_phase = NaN(numel(lfp_ts), n_trials);
choice_phase = NaN(numel(lfp_ts), n_trials);

% bandpass the lfp
f_lfp = bandpass(lfp,freq_range,fs);

% get instantaneous amplitudes and frequencies
[phase,ampenv]= ax_hilbert(f_lfp);

% chop into trials

for t = 1:n_trials
    
    if trials2use(t) == 1
        
        t_cue_time = align_times(t, 1);
        t_choice_time = align_times(t, 2);
        
        cue_erp(:,t) = lfp(t_cue_time - offsets(1) : t_cue_time + offsets(2)-1);
        choice_erp(:,t) = lfp(t_choice_time - offsets(1) : t_choice_time + offsets(2)-1);
        
        cue_amp(:,t) = ampenv(t_cue_time - offsets(1) : t_cue_time + offsets(2)-1);
        choice_amp(:,t) = ampenv(t_choice_time - offsets(1) : t_choice_time + offsets(2)-1);
        
        cue_phase(:,t) = phase(t_cue_time - offsets(1) : t_cue_time + offsets(2)-1);
        choice_phase(:,t) = phase(t_choice_time - offsets(1) : t_choice_time + offsets(2)-1);
        
    end
    
    
end % of looping over trials



end % of function