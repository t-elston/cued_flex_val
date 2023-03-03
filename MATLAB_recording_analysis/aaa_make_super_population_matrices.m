% aaa_make_super_population_matrices.m

% The goal of this function is to make pseudo-population analysis easier.

% This function organizes the firing rates from individual single units into one huge matrix of size n_units x n_time_steps.
% It also generates a corresponding table with the relevant features of each trial. 

datadir = 'C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/recdata/';
savedir = 'C:\Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/population_data/';

rec_files = dir([datadir '*.mat']);

all_FRs = [];
all_OFC_ix = [];
unit_ids = {}; 
trial_features = [];

for f = 1:numel(rec_files)
    
    mat_fname = [datadir rec_files(f).name];
    
    fprintf('\n')  
    disp(['***loading ' rec_files(f).name]) 
    load(mat_fname);

    firing_rates = choice_zFRs;
    ts = choice_down_times;
    
    [n_trials, n_times, n_units] = size(down_choice_FRs);
    
    f_units = [];
    f_OFC_ix = [];
    f_unit_trial_features = [];
    f_unit_ids = {};
    
    % now start constructing the arrays for each unit in the file
    for u = 1:n_units
        
        unit_region=[];
        
        f_units = [f_units ; firing_rates(:,:,u)];
        
        f_unit_trial_features = [f_unit_trial_features ; bhv];
            
        unit_id(1:n_trials,1) = unit_names(u);
        unit_region(1:n_trials,1) = OFC_ix(u);
        
        f_unit_ids = [f_unit_ids ; unit_id];
        
        f_OFC_ix = [f_OFC_ix; unit_region];
       
    end % of looping over units
    
    all_FRs = [all_FRs ; f_units];
    unit_ids = [unit_ids ; f_unit_ids];
    trial_features = [trial_features; f_unit_trial_features];
    all_OFC_ix = [all_OFC_ix; f_OFC_ix];
    
end % of looping over files

% subselect and save the data from each brain area
HPC_data=struct;
OFC_data=struct;

OFC_data.zFRs = all_FRs(all_OFC_ix==1,:);
OFC_data.u_names = unit_ids(all_OFC_ix == 1);
OFC_data.bhv = trial_features(all_OFC_ix==1,:);
OFC_data.bhv_varnames = bhv_varnames;
OFC_data.ts = ts;

HPC_data.zFRs = all_FRs(all_OFC_ix==0,:);
HPC_data.u_names = unit_ids(all_OFC_ix == 0);
HPC_data.bhv = trial_features(all_OFC_ix==0,:);
HPC_data.bhv_varnames = bhv_varnames;
HPC_data.ts = ts;

% now let's save the data
save([savedir, 'OFC_popdata.mat'],'-struct','OFC_data','-v7.3');
save([savedir, 'HPC_popdata.mat'],'-struct','HPC_data','-v7.3');

    





