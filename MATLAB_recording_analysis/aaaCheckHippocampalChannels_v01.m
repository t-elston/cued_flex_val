% aaaCheckHippocampalChannels_v01

datadir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\recording_data\';
LFP_dir = 'E:\CuedFlexVal\sorted_data\';

rec_files = dir([datadir '*.mat']);
LFP_files = dir([LFP_dir '*.pl2']);
LFP_fnames = struct2cell(LFP_files);
LFP_fnames = LFP_fnames(1,:);

n_files = numel(rec_files);



fs = 1000;
[notch_filts] = make_notch_filts(fs);


all_ISIs=[];
all_HPC_ix = [];
u_ctr=0;
good_unit_ctr = 0;
for f = 1:n_files
    
    % find the relevant lfp file name
    lfp_fname_ix = contains(LFP_fnames, rec_files(f).name(1:12));
    lfp_fname = LFP_fnames{lfp_fname_ix};
    

    % behavior
    bhv = load([datadir rec_files(f).name], 'bhv');
    bhv = bhv.bhv;
    
    % sorting notes
    notes = load([datadir rec_files(f).name], 'sorting_notes');
    notes = notes.sorting_notes;
     
    % ISIs
    unit_ISIs = load([datadir rec_files(f).name], 'unit_ISIs');
    unit_ISIs = unit_ISIs.unit_ISIs;
    
    
    
    
    % now let's go through the HPC channels and check for theta and SWRs
    
    HPC_ix = find(contains(notes.Target,'HPC'));
    HPC_notes = notes(HPC_ix,:);
    
    % get details about the LFP channel names
    [~,names] = plx_adchan_names([LFP_dir lfp_fname]);
    names = cellstr(names);
    names = names(cellfun(@(x) ~isempty(strfind(x,'FP')),names)); % only keep FP channels
    
    figure;
    title(rec_files(f).name, 'Interpreter', 'none');
    
    % first, cycle over the units and check their ISIs via the rasters
    for u = 1:numel(HPC_ix)
        u_ctr = u_ctr+1;
        
        u_ISIs = unit_ISIs{HPC_ix(u)};
        
        all_ISIs(u_ctr) = mean(u_ISIs);
        
        
        % now load/ get the LFP trace for the channel this unit was on
        unit_channel = HPC_notes.Channel(u);
        
        % load the LFP for this channel
        [lfp_signal] = load_and_filter_lfp([LFP_dir lfp_fname],names{unit_channel}, notch_filts);
        
        
        fs=1000; % sampling freq
        window=2*fs;
        noverlap = [];
 
        % Calculate a power spectrum with Welch's method
        try
        [psd, freqs] = pwelch(lfp_signal, fs, noverlap, [1 : 50], fs);
        psd = 10*log10(psd);
        
        psd = psd';
        freqs = freqs';
        catch
            xx=[];
        end
        
        ap_mdl = fit(freqs,psd,'exp2');
        
        subtracted_psd = psd - ap_mdl(freqs);
        
        [pks,locs,w,p] = findpeaks(subtracted_psd);
        
        
        if any(locs >= 4 & locs <= 12)
            has_theta = true;
        else
            has_theta = false;
        end
        good_unit_ctr = good_unit_ctr + has_theta;
        
        hold on
        if has_theta
            plot(freqs, zscore(subtracted_psd), 'color', 'r', 'LineWidth',1);
        else
            plot(freqs, zscore(subtracted_psd), 'color', [.5 .5 .5], 'LineWidth',1);
        end

        %ylim([-1,4]);
        
        
    end % of cycling over units
    
    
    

end % of cycling through files

xx=[];

