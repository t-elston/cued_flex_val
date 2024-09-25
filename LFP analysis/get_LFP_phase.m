function [LFP_out, ch_details] = get_LFP_phase(sdata, pl2_fname, pl2_ch_names, pl2_event_times, FOI)

LFP_out = struct; 

% initialize output
LFP_out.theta_phase = [];
LFP_out.theta_amp = [];
LFP_out.beta_phase  = [];
LFP_out.gamma_phase = [];

LFP_out.theta_PLV = [];
LFP_out.beta_PLV = [];
LFP_out.gamma_PLV = [];
LFP_out.alpha_PLV = [];

LFP_out.raw_wave = [];

ch_details = table;
warning('off','MATLAB:table:RowsAddedExistingVars');

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
    
    ch_details.Target(u) = u_probe_notes.Target;
    ch_details.Tower(u) = u_probe_notes.Tower;
    
    % get the lfp channels that were on the same probe as the unit
    candidate_channels = lfp_ch_nums((lfp_ch_nums >= u_probe_notes.First_ch) & (lfp_ch_nums <= u_probe_notes.Last_ch));
    
    % which of those channels on the probe were nearest to the unit?
    [~, closest_ch_ix] = min(abs(candidate_channels - u_channels(u)));
    
    % get the number of that closest channel
    closest_ch = candidate_channels(closest_ch_ix);
    
    % get the index, within the lfp_ch_nums array, that corresponds to that closest channel
    ch_ix = min(find(lfp_ch_nums == closest_ch));
    
    ch_details.ch_name(u) = pl2_ch_names(ch_ix);
    
    % load the LFP
    try
    ChLFP = PL2Ad(pl2_fname, pl2_ch_names{ch_ix});
    [raw_lfp_trials] = chop_lfp_into_trials(ChLFP.Values, pl2_event_times, 2000);
    
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
            
            % now compute the phase alignment for this channel
            [phase_alignment] = compute_LFP_phase_alignment_v2(lfp_trials);
            
            % accumulate the data into fields of a struct
            switch fq
                
                case 1
                    
                    fq_amp = abs(hilbert(adjusted_sig));
                    % smooth the amplitude
                    smoothed_amp = gaussian_smooth_data(fq_amp, 50);
                    [lfp_amp] = chop_lfp_into_trials(smoothed_amp, pl2_event_times, 2000);
          
                    LFP_out.theta_phase = cat(3, LFP_out.theta_phase, lfp_trials);
                    LFP_out.theta_amp = cat(3, LFP_out.theta_amp, lfp_amp);
                    LFP_out.theta_PLV = cat(1, LFP_out.theta_PLV, phase_alignment);
                    
                case 2
                    LFP_out.beta_PLV = cat(1, LFP_out.beta_PLV, phase_alignment);
                case 3
                    LFP_out.gamma_PLV = cat(1, LFP_out.gamma_PLV, phase_alignment);
                case 4
                    LFP_out.alpha_PLV = cat(1, LFP_out.alpha_PLV, phase_alignment);
            end
             
        end % of looping over filters / frequencies

        LFP_out.raw_wave = cat(3, LFP_out.raw_wave, raw_lfp_trials);
        
    end % of looping over units
        
end % of function