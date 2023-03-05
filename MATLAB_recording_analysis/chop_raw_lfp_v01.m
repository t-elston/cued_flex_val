function [out_lfp, lfp_ts] = chop_raw_lfp_v01(lfp, fs, align_times, offsets)

%-------------------------------------------------------------------------
% chops lfp into trials based on the timestamps in align_times and
% trials2use
% Thomas Elston
% 14 Sept 2022
%-------------------------------------------------------------------------
% INPUTS
% lfp         - datastream of unfiltered local field potential
% fs          - sampling frequency of lfp (e.g. 1000 = 1000Hz)
% align_times - n_trials x 1 array of times to align data to 
% offsets     - 2-element array detailing how many milliseconds 
%               [before,after] align_times to use when constructing the LFP scalogram

% OUTPUTS
% out_lfp      - an n_times x n_trials x n_times array aligned to align_times
% lfp_ts       - an n_times array corresponding to the time relative to
%                the aligning event
%-------------------------------------------------------------------------

align_times(align_times==0)=[];

n_trials = numel(align_times);

% compute time stamps
lfp_ts =  -1*offsets(1): offsets(2)-1;

% initialize output
out_lfp = NaN(numel(lfp_ts), n_trials);


% chop into trials
for t = 1:n_trials
        
        t_align_time = align_times(t, 1);
        
        out_lfp(:,t) = lfp(t_align_time - offsets(1) : t_align_time + offsets(2)-1);

end % of looping over trials



end % of function