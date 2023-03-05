% aaa_get_preprocess_LFP_data_main_v01.m

pl2dir = 'E:\CuedFlexVal\sorted_data\';
bhvdir = 'E:\CuedFlexVal\bhv2_files\';
notesdir = 'E:\CuedFlexVal\sorting_notes\';
savedir = 'C:\Users\Thomas Elston\Documents\PYTHON\CuedFlexVal\recdata\';

pl2_files= dir([pl2dir '*.pl2']);
bhv_files = dir([bhvdir '*.bhv2']);
notes_files = dir([notesdir '*.xlsx']);
[n_pl2_files,~] = size(pl2_files);


% some key event codes
keep_event  = 38; % event code indicating a choice was made and the trial should be retained

start_event = 9;
end_event   = 18;
cue_event = 3; % state cues appears on screen
choice_event = 40; % choice options appear on screen
cue_on = 3;
cue_psth_pos_offset = 500;
cue_psth_neg_offset = 500;
choice_psth_pos_offset = 1000;
choice_psth_neg_offset = 1200;


f_ch_tbl = table;
ch_ctr = 0;
all_trial_data={};
all_iti_data  ={};
all_PDC=[];



warning('off','MATLAB:table:RowsAddedExistingVars');
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

for f = 1:n_pl2_files
    
    pl2_fname = [pl2dir pl2_files(f).name];
    
    fprintf('\n')
    disp(['***loading file: ' pl2_files(f).name])
    
    % get indices of trial starts, ends, cue times, and choice times in
    % milliseconds (absolute recording time)
    [trialinfo, Trodalness, Duration] = get_trialindices(pl2_fname, start_event, end_event, cue_event, choice_event, keep_event);
    
    % now I want to make spike tables
    % get some info about the units in this file
    [tscounts, wfcounts, evcounts, slowcounts] = plx_info(pl2_fname,1);
    
    % find the channels with spikes and the number of units on them
    active_channels = find(sum(tscounts) > 0) - 1;
    
    % take trodality into account - if it's stereotrode, only use the odd
    % numbered channels
    if numel(Trodalness) == 2
        channels_to_check = active_channels(logical(mod(active_channels, numel(Trodalness))));
    else
        channels_to_check = active_channels;
    end
    
    n_units_per_channel = sum(tscounts(2:end,channels_to_check+1) > 0);
    total_n_units = sum(n_units_per_channel);
    
    raw_timesteps = -1*cue_psth_neg_offset: cue_psth_pos_offset-1;
    n_trials = numel(trialinfo(:,1));
    
    % load the sorting notes (to get which channel was in which brain area)
    sorting_notes = readtable([notesdir notes_files(f).name]);
    
    % which channels had units on them?
    lfp_channels_to_check = channels_to_check(find(n_units_per_channel > 0));
    
    % extract the LFP and find channels with good signal
    disp('processing LFP...')
    
    % get details about the LFP channel names
    [~,file_ch_names] = plx_adchan_names(pl2_fname);
    file_ch_names = cellstr(file_ch_names);
    file_ch_names = file_ch_names(cellfun(@(x) ~isempty(strfind(x,'FP')),file_ch_names)); % only keep FP channels
    f_ch_tbl= table;
    
    file = cell(numel(lfp_channels_to_check),1);
    brain_area = cell(numel(lfp_channels_to_check),1);
    probe = NaN(numel(lfp_channels_to_check),1);
    ch_num = NaN(numel(lfp_channels_to_check),1);
    ch_name = cell(numel(lfp_channels_to_check),1);
    has_theta = NaN(numel(lfp_channels_to_check),1);
    theta_power = NaN(numel(lfp_channels_to_check),1);


    % find out which channels are HPC and which are OFC
    pw = PoolWaitbar(numel(lfp_channels_to_check), 'computing PSDs...');
    parfor ch_ix = 1:numel(lfp_channels_to_check)
        increment(pw);
        
        % compute PSD and fooof
        ChLFP = PL2Ad(pl2_fname, file_ch_names{lfp_channels_to_check(ch_ix)});
        
        [psd, freqs] = pwelch(ChLFP.Values, 1000, [], [1:50], 1000);
        
        % FOOOF settings
        settings = struct();  % Use defaults
        f_range = [1, 50];
        
        % Run FOOOF, also returning the model
        fooof_results = fooof(freqs', psd', f_range, settings, true);
        new_psd = fooof_results.power_spectrum - fooof_results.ap_fit;
        
        file{ch_ix} = pl2_files(f).name;
        brain_area(ch_ix) = sorting_notes.Target(min(find(sorting_notes.Channel == lfp_channels_to_check(ch_ix))));
        probe(ch_ix) = sorting_notes.TowerNumber(min(find(sorting_notes.Channel == lfp_channels_to_check(ch_ix))));
        ch_num(ch_ix) = lfp_channels_to_check(ch_ix);
        ch_name(ch_ix) = file_ch_names(lfp_channels_to_check(ch_ix));
        has_theta(ch_ix) = any(any(fooof_results.peak_params > 4 & fooof_results.peak_params < 8));
        theta_power(ch_ix) = mean(new_psd(4:8));
        
        
        ChLFP=[];
    end % of looping over channels and finding which have theta
    
    f_ch_tbl.file = file;
    f_ch_tbl.brain_area = brain_area;
    f_ch_tbl.probe = probe;
    f_ch_tbl.ch_num = ch_num;
    f_ch_tbl.ch_name = ch_name;
    f_ch_tbl.has_theta = has_theta;
    f_ch_tbl.theta_power = theta_power;
    
    delete(pw);
    
    
    % only use the channels that have theta peaks
   % f_ch_tbl(f_ch_tbl.has_theta==0,:)=[];
    
    % use only the channel with the biggest theta peak per probe
    probe_ids = unique(f_ch_tbl.probe);
    channels2use = table;
    for p = 1:numel(probe_ids)
        
        this_probe_data = f_ch_tbl(f_ch_tbl.probe == probe_ids(p),:);
        [~,max_theta_ix] = max(this_probe_data.theta_power);
        channels2use(p,:) = this_probe_data(max_theta_ix,:);
        
    end % of looping over probes
    
    
    % loop over channel pairs in this session and compute relevant
    % coherograms
    if ~isempty(channels2use)
        [trial_data, iti_data, to, outPDC] = CFV_compute_coherogram(channels2use, file_ch_names, pl2_fname,  pl2_files(f).name, trialinfo);
        
        all_trial_data = [all_trial_data ; trial_data];
        all_iti_data = [all_iti_data  ; iti_data];
        all_PDC = cat(3,all_PDC,outPDC);
    end
     
end % of looping over pl2 files

% get confidence intervals
for i = 1:numel(all_PDC(1,:,1))
    [PDC_mean(:,i), PDC_sem(:,i)] = GetMeanCI(squeeze(all_PDC(:,i,:))','bootstrap');
end % of cycling over channel pairs



for i = 1:numel(all_trial_data)
[coh_z_map(:,:,i), raw_coh_map(:,:,i)] = CFV_make_coh_map_v01(all_trial_data{i},all_iti_data{i});
end

theta_coh = squeeze(mean(raw_coh_map(4:8,:,:),1))';
beta_coh  = squeeze(mean(raw_coh_map(14:30,:,:),1))';
gamma_coh  = squeeze(mean(raw_coh_map(30:50,:,:),1))';
n_pairs = numel(theta_coh(:,1));

CohFig = figure;
set(CohFig,'renderer','Painters');
set(CohFig,'Units','centimeters','Position',[10 10 40 10]);
subplot(1,3,1);
[t,f] = meshgrid(to,1:50);
coh = surf(t,f,mean(raw_coh_map,3));
% hold on
% plot3([0,0],[1,50],[20,20],'k','LineWidth',2)
% plot3([-700,-700],[1,50],[20,20],'k','LineWidth',2)
view(0,90);
ylabel('Freq (Hz)');
xlabel('Time from Choice');
title('Mean HPC-OFC Coherence')
shading interp
cb = colorbar;
cb.FontSize = 12;
cb.Position = [.16 .66 .01 .2];
%cb.Color = [ 1 1 1 ];
%cb.Box = 'off';
axis tight
set(gca,'FontSize',12,'LineWidth',1);
xlim([-1200,1500])


subplot(1,3,2);
colors = [0 0.5 1; 0.5 0 1; 0.7 0.7 0.7];
hold on
plot(to,mean(theta_coh),'LineWidth',2,'color',colors(1,:));
plot(to,mean(beta_coh),'LineWidth',2,'color',colors(2,:));
plot(to,mean(gamma_coh),'LineWidth',2,'color',colors(3,:));

shadedErrorBar(to,mean(theta_coh),std(theta_coh)/sqrt(n_pairs),'lineprops',{'color',colors(1,:),'LineWidth',2},'patchSaturation',0.4);  
shadedErrorBar(to,mean(beta_coh),std(beta_coh)/sqrt(n_pairs),'lineprops',{'color',colors(2,:),'LineWidth',2},'patchSaturation',0.4);  
shadedErrorBar(to,mean(gamma_coh),std(gamma_coh)/sqrt(n_pairs),'lineprops',{'color',colors(3,:),'LineWidth',2},'patchSaturation',0.4);
set(gca,'FontSize',12,'LineWidth',1);
xlim([-1200,1500])
ylabel('Coherence');
xlabel('Time from Pics On');


CT = cbrewer('qual','Set1',9);
subplot(1,3,3);
plot([1:50],PDC_mean(:,5),'color',CT(2,:),'LineWidth',3);
plot([1:50],PDC_mean(:,6),'color',[.1,.1,.1],'LineWidth',3);

shadedErrorBar([1:50],PDC_mean(:,5), PDC_sem(:,5),'lineprops',{'color',CT(2,:),'LineWidth',3},'patchSaturation',0.4);
shadedErrorBar([1:50],PDC_mean(:,6), PDC_sem(:,6),'lineprops',{'color',[.1,.1,.1],'LineWidth',3},'patchSaturation',0.4);
ylim([-.04, .15]);
xlabel('Frequency (Hz)');
ylabel('Net Signal Directionality');
legend('real signal','shuffled signal');
set(gca,'FontSize',12,'LineWidth',1);







