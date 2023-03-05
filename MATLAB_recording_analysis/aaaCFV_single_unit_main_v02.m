% aaaCFV_single_unit_main_v01

% where are the data located?
datadir = 'C:\Users\Thomas Elston\Documents\PYTHON\CuedFlexVal\recdata\';

rec_files = dir([datadir '*.mat']);

n_files = numel(rec_files);

% initialize some arrays for aggregating
all_sig_factors=[];
all_R2 = [];
all_ofc_ix =[];

for f = 1:n_files
    
    % load the data for this file
    fprintf('\n');
    fprintf(['***assessing file: ' rec_files(f).name ' \n']);
    
    % behavior
    bhv = load([datadir rec_files(f).name], 'bhv'); bhv = array2table(bhv.bhv);
    bhv_varnames = load([datadir rec_files(f).name], 'bhv_varnames'); bhv_varnames=bhv_varnames.bhv_varnames;
    bhv.Properties.VariableNames = bhv_varnames;
    
    % firing_rates
    FRs = load([datadir rec_files(f).name], 'down_choice_FRs');
    FRs = FRs.down_choice_FRs;

    % rasters and timesteps
    rasters = load([datadir rec_files(f).name], 'choice_rasters');
    ts = load([datadir rec_files(f).name], 'choice_down_times');
    ts = ts.choice_down_times;
    rasters = rasters.choice_rasters;
     
    % brain area indices
    ofc_ix = load([datadir rec_files(f).name], 'OFC_ix');
    ofc_ix = ofc_ix.OFC_ix;
    
    OFC_ix = ofc_ix==1;
    HPC_ix = ofc_ix==0;
    
    [n_trials, n_times, n_units] = size(FRs);
    
    sig_factors = [];
    shuff_sig_factors=[];
    R2=[];
    factors = {bhv.state, bhv.chosenval, bhv.side};
    fnames = {'state','val','side'};
    % now loop over units and timesteps
    pw = PoolWaitbar(n_units, 'running selectivity analysis...');
    parfor u = 1:n_units
        increment(pw);
        
        %[sig_factors(:,:,u), shuff_sig_factors(:,:,u)] = parallel_sliding_anova(FRs(:,:,u), factors,fnames);  
        [sig_factors(:,:,u), R2(:,:,u)] = parallel_sliding_GLM(FRs(:,:,u), bhv, .05); 
        % make a PSTH and raster plot
%         if cue_sig(u_ctr,1) == 1
%             PSTH_raster_plot_v02(bhv,FRs(:,:,u),ts,rasters(:,:,u), ts, ofc_ix(u,:),cue_pvals(u_ctr,:), choice_pvals(u_ctr,:));
%         end
    end % of looping over units
    delete(pw);
     
    all_sig_factors = cat(3,all_sig_factors , sig_factors);
    all_R2 = cat(3,all_R2,R2);
    all_ofc_ix = [all_ofc_ix;OFC_ix];
    
end % of looping over files

state_sig = squeeze(all_sig_factors(1,:,:))';
val_sig = squeeze(all_sig_factors(2,:,:))';

figure;
subplot(1,2,1);
hold on
plot(ts,mean(state_sig(all_ofc_ix ==0,:)))
plot(ts,mean(val_sig(all_ofc_ix ==0,:)))

subplot(1,2,2);
hold on
plot(ts,mean(state_sig(all_ofc_ix ==1,:)))
plot(ts,mean(val_sig(all_ofc_ix ==1,:)))





xx=[];