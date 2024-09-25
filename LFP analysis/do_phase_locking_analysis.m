function [cue_pvals, cue_PLV, cue_mean_angle, ...
          pics_pvals, pics_PLV, pics_mean_angle] = do_phase_locking_analysis(sdata, pl2_fname, pl2_ch_names, pl2_event_times, FOI)

% construct a cell array containing all of the bandpass filters
filters = {};

for f = 1:numel(FOI(:,1))
    
    l_bound = FOI(f,1);
    
    u_bound = FOI(f,2);
     
    filters{f} = designfilt('bandpassfir', ... response type
    'FilterOrder',1000, ... filter order
    'CutoffFrequency1', l_bound,... lower range
    'CutoffFrequency2', u_bound,... upper range
    'SampleRate', 1000); ... sampling rate
       
end % of looping over filters to construct

% keep only the numeric part of the pl2_ch_names
for lfp_ch = 1:numel(pl2_ch_names)
    
    lfp_ch_nums(lfp_ch,1) = str2num(pl2_ch_names{lfp_ch}(end-2 :end));
    
end % of looping over lfp channel names

% find the relevant channels to assess
u_channels = sdata.unit_ch_nums;

% loop over each unit's channel
for u = 1:numel(u_channels)
    
    % find the nearest relevant LFP channel
    u_probe = sdata.unit_probe_nums(u);
    u_probe_notes = sdata.probe_notes(sdata.probe_notes.Tower == u_probe, :);
    
    % get the lfp channels that were on the same probe as the unit
    candidate_channels = lfp_ch_nums((lfp_ch_nums >= u_probe_notes.First_ch) & (lfp_ch_nums <= u_probe_notes.Last_ch));
    
    % which of those channels on the probe were nearest to the unit?
    [~, closest_ch_ix] = min(abs(candidate_channels - u_channels(u)));
    
    % get the number of that closest channel
    closest_ch = candidate_channels(closest_ch_ix);
    
    % get the index, within the lfp_ch_nums array, that corresponds to that closest channel
    ch_ix = min(find(lfp_ch_nums == closest_ch));
    
    % load the LFP
    try
    ChLFP = PL2Ad(pl2_fname, pl2_ch_names{ch_ix});
    catch
        xx=[];
    end
    

        % loop over the different bp filters
        for fq = 1:numel(filters)
            
            fq_filt = filters{fq};
            
            % apply the bandpass
            fq_signal = fftfilt(fq_filt.Coefficients,ChLFP.Values);
            
            % adjust time to account for filter lag
            LFPtime = 1:numel(ChLFP.Values); %ms
            lag = 1+(1+1)*ChLFP.ADFreq/2;
            adjusted_sig =fq_signal(lag:end);
            
            % extract the instantaneous phase
            fq_phase = angle(hilbert(adjusted_sig));
            
            % now chop the LFP phase into trials
            [lfp_trials] = chop_lfp_into_trials(fq_phase, pl2_event_times, 2000);
            
            [~, cue_start] = min(abs(sdata.raster_ts - -700));
            [~, pics_start] = min(abs(sdata.raster_ts - 0));
            
            
            % look at clustering around the circle during the cue period
            [cue_pvals(u, fq), cue_PLV(u, fq), cue_mean_angle(u, fq)]  = assess_circular_clustering(sdata.rasters(:, cue_start: cue_start+500, u),...
                                                                                       lfp_trials(:, cue_start: cue_start+500));
                                                                                   
             % look at clustering around the circle during the cue and choice period
            [pics_pvals(u, fq), pics_PLV(u, fq), pics_mean_angle(u, fq)]  = assess_circular_clustering(sdata.rasters(:, pics_start: pics_start+500, u),...
                                                                                       lfp_trials(:, pics_start: pics_start+500));
            
        end % of looping over filters / frequencies
    end % of looping over units

end % of function