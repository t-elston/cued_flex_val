function [pval, PLV, mean_angle] = assess_circular_clustering(rasters, lfp_phase)

% checks whether unit spikes occur at a consisent phase angle of the lfp
% inputs: 
% spikes = n_trials x n_times x n_units array of rasters 
% lfp_phase = n_trials x n_times x n_channels of instaneous phase data
% of band-passed LFP data derived from the Hilbert transform
% unit_channels = n_units x 1 array of channel numbers; indicates which LFP
% channel the unit is paired with
% lfp_ch_num = n_lfp_channels x 1 array indicating the absolute channel
% number a set of phase data corresponds to


[n_trials, n_times, n_units] = size(rasters);

rasters(isnan(rasters)) = 0;

unit_spikes = logical(rasters);

% get parameters from circ_stats toolbox
pval = circ_rtest(lfp_phase(unit_spikes));
mean_angle = circ_mean(lfp_phase(unit_spikes));
PLV = circ_r(lfp_phase(unit_spikes));



end % of function