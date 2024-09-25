% aaa_make_super_population_matrices_v02.m

% The goal of this function is to make pseudo-population analysis easier.

% This function organizes the firing rates from individual single units into one huge matrix of size n_units x n_time_steps.
% It also generates a corresponding table with the relevant features of each trial. 

datadir = 'C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/recdata2/';
savedir = 'C:\Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/population_data/';

rec_files = dir([datadir '*.mat']);

% subselect and save the data from each brain area
OFC_data=struct;
OFC_data.zFRs=[];
OFC_data.u_names=[];
OFC_data.bhv=[];
OFC_data.n_trials_per_choice_cond=[];


HPC_data=struct;
HPC_data.zFRs=[];
HPC_data.u_names=[];
HPC_data.bhv=[];
HPC_data.n_trials_per_choice_cond=[];


for f = 1:numel(rec_files)
    
    mat_fname = [datadir rec_files(f).name];
    
    fprintf('\n')  
    disp(['***loading ' rec_files(f).name]) 
    load(mat_fname);
    
    % make table of the behavior and count number of trials per condition
    bhv2 = array2table(sdata.bhv,'VariableNames',sdata.bhv_varnames);

    firing_rates = sdata.nrml_FRs;
    t_mids = sdata.t_mids;
    
    [n_trials, n_times, n_units] = size(firing_rates);
    
    % now start constructing the arrays for each unit in the file
    for u = 1:n_units
        
        % was this a HPC or OFC unit?
        if sdata.OFC_ix(u) == 1
            
            OFC_data.zFRs = [OFC_data.zFRs ; firing_rates(:,:,u)];
            OFC_data.u_names = [OFC_data.u_names ; repmat(sdata.unit_names(u),n_trials,1)]; 
            OFC_data.bhv = [OFC_data.bhv ; sdata.bhv];
            
            % ensure there are enough trials/condition during the choice for this neuron
            OFC_data.n_trials_per_choice_cond= cat(2, OFC_data.n_trials_per_choice_cond,... 
                              [ sum(bhv2.state == 1 & bhv2.chosenval == 1) ; 
                                sum(bhv2.state == 1 & bhv2.chosenval == 2) ;
                                sum(bhv2.state == 1 & bhv2.chosenval == 3) ;
                                sum(bhv2.state == 1 & bhv2.chosenval == 4) ;
                                sum(bhv2.state == -1 & bhv2.chosenval == 1) ; 
                                sum(bhv2.state == -1 & bhv2.chosenval == 2) ;
                                sum(bhv2.state == -1 & bhv2.chosenval == 3) ;
                                sum(bhv2.state == -1 & bhv2.chosenval == 4)]) ;
            
        else
            
            HPC_data.zFRs = [HPC_data.zFRs ; firing_rates(:,:,u)];
            HPC_data.u_names = [HPC_data.u_names ; repmat(sdata.unit_names(u),n_trials,1)]; 
            HPC_data.bhv = [HPC_data.bhv ; sdata.bhv];
            
            % ensure there are enough trials/condition during the choice for this neuron
            HPC_data.n_trials_per_choice_cond= cat(2, HPC_data.n_trials_per_choice_cond,... 
                              [ sum(bhv2.state == 1 & bhv2.chosenval == 1) ; 
                                sum(bhv2.state == 1 & bhv2.chosenval == 2) ;
                                sum(bhv2.state == 1 & bhv2.chosenval == 3) ;
                                sum(bhv2.state == 1 & bhv2.chosenval == 4) ;
                                sum(bhv2.state == -1 & bhv2.chosenval == 1) ; 
                                sum(bhv2.state == -1 & bhv2.chosenval == 2) ;
                                sum(bhv2.state == -1 & bhv2.chosenval == 3) ;
                                sum(bhv2.state == -1 & bhv2.chosenval == 4)]) ;
                            
        end % of determining whether this is a HPC or OFC unit
               
    end % of looping over units
    
end % of looping over files

OFC_data.bhv_varnames = sdata.bhv_varnames;
OFC_data.ts = t_mids;

HPC_data.bhv_varnames = sdata.bhv_varnames;
HPC_data.ts = t_mids;

% now let's save the data
save([savedir, 'OFC_popdata.mat'],'-struct','OFC_data','-v7.3');
save([savedir, 'HPC_popdata.mat'],'-struct','HPC_data','-v7.3');

    





