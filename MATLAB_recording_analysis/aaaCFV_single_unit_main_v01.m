% aaaCFV_single_unit_main_v01

% where are the data located?
datadir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\CuedFlexVal\recording_data\';

rec_files = dir([datadir '*.mat']);

n_files = numel(rec_files);

% initialize some arrays for aggregating
all_sig_factors=[];
all_s_sig_factors=[];
all_pE2=[];
all_s_pE2 = [];
all_OFC_ix =[];
all_HPC_ix = [];
u_ctr = 0;
cue_pvals=[];
cue_betas = [];
choice_pvals=[];
choice_betas=[];
for f = 1:n_files
    
    % load the data for this file
    fprintf('\n');
    fprintf(['***assessing file: ' rec_files(f).name ' \n']);
    
    %---
    % behavior
    %---
    bhv = load([datadir rec_files(f).name], 'bhv');
    bhv = bhv.bhv;
    
    %---
    % firing_rates
    %---
    firing_rates = load([datadir rec_files(f).name], 'down_firing_rates');
    firing_rates = firing_rates.down_firing_rates;
    
    
    %---
    % rasters
    %---
    rasters = load([datadir rec_files(f).name], 'rasters');
    raster_ts = load([datadir rec_files(f).name], 'raw_times');
    raster_ts = raster_ts.raw_times;
    rasters = rasters.rasters;
     
    %---
    % timesteps
    %---
    ts = load([datadir rec_files(f).name], 'down_timesteps');
    ts = ts.down_timesteps;
    
    %---
    % sorting notes
    %---
    notes = load([datadir rec_files(f).name], 'sorting_notes');
    notes = notes.sorting_notes;
    
    
    
    OFC_ix = contains(notes.Target, 'OFC');
    HPC_ix = contains(notes.Target, 'HPC');
    
    % only use the correct trials 
    hit_ix = bhv.pickedbest == 1; 
    hilo_ix = bhv.chosenval < 3;
    firing_rates = firing_rates(bhv.UseTrial == 1,:,:);
    bhv = bhv(bhv.UseTrial ==1,:);
    
    % create the factors for the anova 
    state = bhv.state;
    val = bhv.chosenval;

    forced = bhv.forcedchoice ==1; 
    hilo_ix = bhv.chosenval < 3;
    type = bhv.state_type;
    state_val = state.*val;
    side = bhv.side;
    picked_best = bhv.pickedbest;
    chosen_im = bhv.chosen_im;
    %factors = {state, val, chosen_im};
    factors = {state, val, side};
    fnames = {'state','val','side'};
    
    [n_trials, n_times, n_units] = size(firing_rates);
    
    sig_factors = [];
    shuff_sig_factors=[];
    % now loop over units and timesteps
    pw = PoolWaitbar(n_units, 'running selectivity analysis...');
    for u = 1:n_units
         increment(pw);
         u_ctr = u_ctr+1;

            %[sig_factors(:,:,u), shuff_sig_factors(:,:,u)] = parallel_sliding_GLM(firing_rates(:,:,u), bhv, .05);
        % make a PSTH and raster plot
        %PSTH_raster_plot(bhv,firing_rates(:,:,u),ts,rasters(:,:,u), raster_ts, notes(u,:),sig_factors(:,:,u));

    end % of looping over units
     delete(pw);

     
     % temporally threshold selectivity
     %[new_sig_factors] = temporally_threshold_selectivity(sig_factors, 8);
     %[new_sig_s_factors] = temporally_threshold_selectivity(shuff_sig_factors, 2);

    
%      % aggregate this file's data into a larger array
%      all_sig_factors = cat(3, all_sig_factors, sig_factors);
%      all_s_sig_factors = cat(3, all_s_sig_factors, shuff_sig_factors);
%      
    
    all_OFC_ix = [all_OFC_ix ; OFC_ix];
    all_HPC_ix = [all_HPC_ix ; HPC_ix];
    
end % of looping over files

% val_sig = squeeze(all_sig_factors(1,:,:));
% 
% cmap = cbrewer('qual','Paired',12);
% figure;
% subplot(1,2,1);
% hold on
% plot(ts,nanmean(all_sig_factors(1,:,logical(all_OFC_ix)),3), 'LineWidth',2, 'color',cmap(2,:));
% plot(ts,nanmean(all_sig_factors(2,:,logical(all_OFC_ix)),3),'LineWidth',2, 'color',cmap(4,:));
% plot(ts,nanmean(all_sig_factors(3,:,logical(all_OFC_ix)),3),'LineWidth',2, 'color',cmap(6,:));
% 
% plot(ts,nanmean(all_s_sig_factors(1,:,logical(all_OFC_ix)),3), 'LineWidth',2, 'color',cmap(1,:));
% plot(ts,nanmean(all_s_sig_factors(2,:,logical(all_OFC_ix)),3),'LineWidth',2, 'color',cmap(3,:));
% plot(ts,nanmean(all_s_sig_factors(3,:,logical(all_OFC_ix)),3),'LineWidth',2, 'color',cmap(5,:));
% 
% plot([0 0], ylim,'k','LineWidth',2);
% plot([-700 -700], ylim,'k','LineWidth',2);
% xlim([-1000 1500]);
% 
% 
% legend('State','Value','State*Value');
% ylabel('prop. sig units');
% xlabel('time from choice on');
% title('OFC');
% 
% subplot(1,2,2);
% hold on
% plot(ts,nanmean(all_sig_factors(1,:,logical(all_HPC_ix)),3), 'LineWidth',2, 'color',cmap(2,:));
% plot(ts,nanmean(all_sig_factors(2,:,logical(all_HPC_ix)),3),'LineWidth',2, 'color',cmap(4,:));
% plot(ts,nanmean(all_sig_factors(3,:,logical(all_HPC_ix)),3),'LineWidth',2, 'color',cmap(6,:));
% 
% plot(ts,nanmean(all_s_sig_factors(1,:,logical(all_HPC_ix)),3), 'LineWidth',2, 'color',cmap(1,:));
% plot(ts,nanmean(all_s_sig_factors(2,:,logical(all_HPC_ix)),3),'LineWidth',2, 'color',cmap(3,:));
% plot(ts,nanmean(all_s_sig_factors(3,:,logical(all_HPC_ix)),3),'LineWidth',2, 'color',cmap(5,:));
% title('HPC');
% 
% plot([0 0], ylim,'k','LineWidth',2);
% plot([-700 -700], ylim,'k','LineWidth',2);
% xlim([-1000 1500]);
% 
% 





xx=[];