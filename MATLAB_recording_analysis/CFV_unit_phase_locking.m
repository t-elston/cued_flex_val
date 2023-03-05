function [phase_locking_results] = CFV_unit_phase_locking(spikes, lfp_phase, lfp_wave, unit_channels, lfp_ch_num)

% checks whether unit spikes occur at a consisent phase angle of the lfp
% inputs: 
% spikes = n_trials x n_times x n_units array of rasters 
% lfp_phase = n_trials x n_times x n_channels of instaneous phase data
% of band-passed LFP data derived from the Hilbert transform
% unit_channels = n_units x 1 array of channel numbers; indicates which LFP
% channel the unit is paired with
% lfp_ch_num = n_lfp_channels x 1 array indicating the absolute channel
% number a set of phase data corresponds to



[n_trials, n_times, n_units] = size(spikes);

phase_locking_results = NaN(n_units,4);


for u = 1:n_units
    
    unit_lfp_phase = lfp_phase(:,:, find(unit_channels(u) == lfp_ch_num));
    unit_lfp_wave = lfp_wave(:,:, find(unit_channels(u) == lfp_ch_num));

    unit_spikes = logical(spikes(:,:,u));
    
    spike_phases = rad2deg(unit_lfp_phase(unit_spikes));
    [N,edges,bin] = histcounts(spike_phases,18);
    bin_centers = edges(2:end) - .5*mean(diff(edges));
    
    bin_ids = unique(bin);
    bin_phases = NaN(size(bin));
    for b = 1:numel(bin_ids)
        b_ix = bin == bin_ids(b);
    bin_phases(b_ix) = bin_centers(b);
    end
    
    % convert back to radians
    bin_phases = rad2deg(bin_phases);
    spike_phases = rad2deg(spike_phases);

    
    % get parameters from circ_stats toolbox
    [phase_locking_results(u,1), phase_locking_results(u,4)] = circ_rtest(bin_phases);
    phase_locking_results(u,2) = circ_mean(bin_phases);
    phase_locking_results(u,3) = circ_r(bin_phases);    
    
    
end % of looping over units

end % of function