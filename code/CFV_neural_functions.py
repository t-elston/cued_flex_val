'''
CFV_neural_functions.py
Module containing functions called by CFV_neural_main.py
'''
# imports
import glob
from secrets import choice
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import statsmodels.api as sm
import pandas as pd
import pingouin as pg
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import metrics
from statsmodels.distributions.empirical_distribution import ECDF
import numpy.matlib
from scipy import stats
import seaborn as sns
import utils as ut

import decoding_utilities as dut
import mat73

def unit_selectivity_analysis(datadir):
    
    '''
    loop over each neuron in each file and assess its coding profile
    '''

    # initialize some output dataframes
    cue_ttests = pd.DataFrame()
    choice_glms = pd.DataFrame()

    # find file names
    datafiles = glob.glob(datadir + '*.mat')
    all_bhv = pd.DataFrame()

    for fnum, fname in enumerate(datafiles):
        print('file ' + str(fnum+1) + ' \ ' + str(np.size(datafiles)))

        # load this file
        f_data = mat73.loadmat(fname)

        # pull out this file's data
        bhv= pd.DataFrame(f_data['bhv'], columns=f_data['bhv_varnames'])
        cue_phaselocking= pd.DataFrame(f_data['cue_phase_locking'], columns=f_data['phase_locking_varnames'])
        choice_phaselocking= pd.DataFrame(f_data['choice_phase_locking'], columns=f_data['phase_locking_varnames'])

        cue_spikes = f_data['raw_cue_FRs']
        choice_spikes = f_data['choice_zFRs']
        choice_ts = f_data['choice_down_times']
        cue_ts = f_data['cue_raw_times']
        OFC_ix = f_data['OFC_ix'] == 1
        unit_names = f_data['unit_names']

        # delete the big dict to save some space
        del f_data

        # accumulate all behavior
        all_bhv = all_bhv.append(bhv,ignore_index=True)

        # get the mean firing rates for the GLM
        cue_win_start = np.argmin(np.abs(cue_ts - 10))
        cue_win_end = np.argmin(np.abs(cue_ts - 400)) 
        choice_win_start = np.argmin(np.abs(choice_ts - 10))
        choice_win_end = np.argmin(np.abs(choice_ts - 400)) 
        fix_win_start = np.argmin(np.abs(cue_ts - -300))
        fix_win_end = np.argmin(np.abs(cue_ts - 10)) 

        fix_mean_frs = np.mean(cue_spikes[:,fix_win_start:fix_win_end,:], axis=1)
        cue_mean_frs = np.mean(cue_spikes[:,cue_win_start:cue_win_end,:], axis=1)
        choice_mean_frs = np.mean(choice_spikes[:,choice_win_start:choice_win_end,:], axis=1)
        n_trials, n_units = cue_mean_frs.shape
        cue_ids = np.zeros(shape=(len(bhv['state'],)))
        cue_ids[(bhv['state']==1) & (bhv['state_type']==1)] =1
        cue_ids[(bhv['state']==1) & (bhv['state_type']==-1)] =2
        cue_ids[(bhv['state']==-1) & (bhv['state_type']==1)] =3
        cue_ids[(bhv['state']==-1) & (bhv['state_type']==-1)] =4

        # get the data into a table
        glm_df = pd.DataFrame()
        glm_df['state'] = bhv['state']
        glm_df['state_type'] = bhv['state_type']
        glm_df['cue_id'] = cue_ids.astype(int)
        glm_df['t_num'] = np.arange(n_trials)
        glm_df['val'] = bhv['chosenval']
        glm_df['side'] = bhv['side']

        # get some useful indices
        t1_ix = glm_df['state_type'] == 1
        t2_ix = glm_df['state_type'] == -1
        s1_ix = glm_df['state'] == 1
        s2_ix = glm_df['state'] == -1

        f_cue_ttest = pd.DataFrame()
        f_choice_glm = pd.DataFrame()

        for u in range(n_units):

            # analyze the choice spikes (z-scored firing rates so betas are comparable)
            glm_df['choice_fr'] = choice_mean_frs[:,u]
            glm_df['cue_fr'] = cue_mean_frs[:,u]
            glm_df['fix_fr'] = fix_mean_frs[:,u]

            choice_mdl = smf.glm(formula = 'choice_fr ~ state*val + t_num',
                            data = glm_df).fit()

            f_choice_glm.at[u,'file']   = fname
            f_choice_glm.at[u,'unit_num'] = u
            f_choice_glm.at[u,'u_name'] = unit_names[u]
            f_choice_glm.at[u,'ofc']   = OFC_ix[u].astype(int)
            f_choice_glm.at[u,'state_p'] = choice_mdl.pvalues['state']
            f_choice_glm.at[u,'val_p'] = choice_mdl.pvalues['val']
            f_choice_glm.at[u,'state:val_p'] = choice_mdl.pvalues['state:val']
            f_choice_glm.at[u,'state_b'] = choice_mdl.params['state']
            f_choice_glm.at[u,'val_b'] = choice_mdl.params['val']
            f_choice_glm.at[u,'state:val_b'] = choice_mdl.params['state:val']
            f_choice_glm.at[u,'rayleigh_p'] = choice_phaselocking['pval'].loc[u]
            f_choice_glm.at[u,'angle'] = choice_phaselocking['mean_angle'].loc[u]
            f_choice_glm.at[u,'vlen'] = choice_phaselocking['vector_len'].loc[u]

            
            # now look at value nested within state
            s1_choice_mdl = smf.glm(formula = 'choice_fr ~ val + t_num',
                            data = glm_df.loc[s1_ix,:]).fit()

            s2_choice_mdl = smf.glm(formula = 'choice_fr ~ val + t_num',
                            data = glm_df.loc[s2_ix,:]).fit()

            f_choice_glm.at[u,'s1_val_p'] = s1_choice_mdl.pvalues['val']
            f_choice_glm.at[u,'s1_val_int'] = s1_choice_mdl.params['Intercept']
            f_choice_glm.at[u,'s1_val_b'] = s1_choice_mdl.params['val']

            f_choice_glm.at[u,'s2_val_p'] = s2_choice_mdl.pvalues['val']
            f_choice_glm.at[u,'s2_val_int'] = s2_choice_mdl.params['Intercept']
            f_choice_glm.at[u,'s2_val_b'] = s2_choice_mdl.params['val']


            # how discriminable is an arbitrary combination?
            s1_t1_frs = glm_df['cue_fr'].loc[s1_ix & t1_ix]
            s1_t2_frs = glm_df['cue_fr'].loc[s1_ix & t2_ix]
            s2_t1_frs = glm_df['cue_fr'].loc[s2_ix & t1_ix]
            s2_t2_frs = glm_df['cue_fr'].loc[s2_ix & t2_ix]

            state_frs = np.concatenate((s1_t1_frs, s1_t2_frs, s2_t1_frs, s2_t2_frs), axis =0)

            state_labels = np.concatenate((np.zeros(len(s1_t1_frs), ),
                                           np.zeros(len(s1_t2_frs), ),
                                           np.ones(len(s2_t1_frs), ),
                                           np.ones(len(s2_t2_frs), )), axis = 0).astype(int)

            type_labels = np.concatenate((np.zeros(len(s1_t1_frs), ),
                                           np.ones(len(s1_t2_frs), ),
                                           np.zeros(len(s2_t1_frs), ),
                                           np.ones(len(s2_t2_frs), )), axis = 0).astype(int)

            arbitrary_labels = np.concatenate((np.zeros(len(s1_t1_frs), ),
                                               np.ones(len(s1_t2_frs), ),
                                               np.ones(len(s2_t1_frs), ),
                                               np.zeros(len(s2_t2_frs), )), axis = 0).astype(int)

            shuffled_labels = np.random.permutation(state_labels).astype(int)

            f_cue_ttest.at[u,'file']   = fname
            f_cue_ttest.at[u,'unit_num'] = u
            f_cue_ttest.at[u,'u_name'] = unit_names[u]
            f_cue_ttest.at[u,'ofc']   = OFC_ix[u].astype(int)
            f_cue_ttest.at[u,'state_p'] = pg.ttest(cue_mean_frs[s1_ix,u],cue_mean_frs[s2_ix,u])['p-val'].values
            f_cue_ttest.at[u,'type_p'] = pg.ttest(state_frs[type_labels==1],state_frs[type_labels==0])['p-val'].values
            f_cue_ttest.at[u,'mix_p'] = pg.ttest(state_frs[arbitrary_labels==1],state_frs[arbitrary_labels==0])['p-val'].values
            f_cue_ttest.at[u,'shuffle_p'] = pg.ttest(state_frs[shuffled_labels==1],state_frs[shuffled_labels==0])['p-val'].values
            f_cue_ttest.at[u,'state_d'] = pg.ttest(state_frs[state_labels==1],state_frs[state_labels==0])['cohen-d'].values
            f_cue_ttest.at[u,'type_d'] =pg.ttest(state_frs[type_labels==1],state_frs[type_labels==0])['cohen-d'].values
            f_cue_ttest.at[u,'mix_d'] = pg.ttest(state_frs[arbitrary_labels==1],state_frs[arbitrary_labels==0])['cohen-d'].values
            f_cue_ttest.at[u,'shuffle_d'] = pg.ttest(state_frs[shuffled_labels==1],state_frs[shuffled_labels==0])['cohen-d'].values

            # do a 1-way anova on the 4 state cues in the fixation period
            f_cue_ttest.at[u,'fix_p'] = pg.anova(data=glm_df,dv='fix_fr',between='cue_id')['p-unc'].values[0]

            f_cue_ttest.at[u,'rayleigh_p'] = cue_phaselocking['pval'].loc[u]
            f_cue_ttest.at[u,'angle'] = cue_phaselocking['mean_angle'].loc[u]
            f_cue_ttest.at[u,'vlen'] = cue_phaselocking['vector_len'].loc[u]


        # aggregate these values into some output dataframes
        cue_ttests = cue_ttests.append(f_cue_ttest,ignore_index=True)
        choice_glms = choice_glms.append(f_choice_glm,ignore_index=True)

    choice_glms['cell_id'] = np.arange(len(choice_glms))
    cue_ttests['cell_id'] = np.arange(len(cue_ttests))

    return cue_ttests, choice_glms, all_bhv


def plot_cue_selective_units(cue_ttests):

    savedir = 'C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/cue_units/'
    # plot rasters and psth aligned to cue on
    state_coding_ix = (cue_ttests['state_p'] < .05)
    state_cells = cue_ttests.loc[state_coding_ix,:]
    state_cells = state_cells.reset_index(drop=True)
    ofc_ix = state_cells['ofc']==1

    # loop over individual cells and plot
    for u in range(len(state_cells)):

        print(str(u) +  ' / ' + str(len(state_cells)))

        # load this unit's file and extract this unit
        f_data = mat73.loadmat(state_cells['file'][u])
        spikes = f_data['raw_cue_FRs']
        rasters = f_data['cue_rasters']
        u_name = state_cells['u_name'][u][0]
        u_in_ofc = state_cells['ofc'][u]==1
        ts = f_data['cue_raw_times']
        ts = ts[0:len(rasters[0,:,0])]
        bhv= pd.DataFrame(f_data['bhv'], columns=f_data['bhv_varnames'])
        del f_data

        # pull out this unit's spikes and rasters
        u_spikes = spikes[:,:,(state_cells['unit_num'][u]).astype(int)]
        u_rasters = rasters[:,:,(state_cells['unit_num'][u]).astype(int)]

        # get indices for the states and cues
        s1_ix = (bhv['state'] ==1)
        s2_ix = (bhv['state'] ==-1)
        t1_ix = (bhv['state_type'] ==1)
        t2_ix = (bhv['state_type'] ==-1)

        # get state-type combination spikes
        s1_t1_spikes = np.mean(u_spikes[s1_ix & t1_ix,:],axis =0)
        s1_t2_spikes = np.mean(u_spikes[s1_ix & t2_ix,:],axis =0)
        s2_t1_spikes = np.mean(u_spikes[s2_ix & t1_ix,:],axis =0)
        s2_t2_spikes = np.mean(u_spikes[s2_ix & t2_ix,:],axis =0)

        # n_trials to show rasters for each condition
        n_trials = 250
        # grab a random n_trials for raster plotting
        s1_t1_rasters = u_rasters[np.random.permutation(np.where(s1_ix & t1_ix)[0])[0:n_trials],:]
        s1_t2_rasters = u_rasters[np.random.permutation(np.where(s1_ix & t2_ix)[0])[0:n_trials],:]
        s2_t1_rasters = u_rasters[np.random.permutation(np.where(s2_ix & t1_ix)[0])[0:n_trials],:]
        s2_t2_rasters = u_rasters[np.random.permutation(np.where(s2_ix & t2_ix)[0])[0:n_trials],:]

        # define a colormap
        cmap = plt.cm.Paired(np.linspace(0,1,12))
        # define the colors for the states/types
        s1_t1_c = cmap[1,:] # sat blue
        s1_t2_c = cmap[0,:] # unsat blue
        s2_t1_c = cmap[5,:] # sat red
        s2_t2_c = cmap[4,:] # unsat red

        # put the rasters into a big array and make a corresponding array of colors
        all_rasters = np.concatenate((s1_t1_rasters,s1_t2_rasters,s2_t1_rasters,s2_t2_rasters), axis=0)
        all_colors  = np.concatenate((np.matlib.repmat(s1_t1_c,n_trials,1),
                                      np.matlib.repmat(s1_t2_c,n_trials,1),
                                      np.matlib.repmat(s2_t1_c,n_trials,1),
                                      np.matlib.repmat(s2_t2_c,n_trials,1)), axis = 0)

        # update the row values of the rasters so we can plot row by row  
        for t in range(len(all_rasters)): 
            all_rasters[t,:] = (all_rasters[t,:]*t) 
        all_rasters[all_rasters == 0] = np.nan   

        # get details of unit name
        if u_in_ofc:
            u_details = 'OFC_' + u_name
        else:
            u_details = 'HPC_' + u_name


        #---------------------------------------
        #               PLOT
        #---------------------------------------
        sns.set_style("ticks")
        fig = plt.figure(figsize=(6, 6), dpi=300)
        grid = fig.add_gridspec(100, 6)
        # plot the rasters
        raster_ax = fig.add_subplot(grid[0:70, 0:6])
        for t in range(len(all_rasters)): 
            plt.scatter(ts,all_rasters[t,:], color = all_colors[t,:], s=1, marker='.')
        plt.gca().axis('off')
        plt.xlim([-300, 500])
        #plt.plot([0,0], plt.gca().get_ylim(), color = 'black')

        psth_ax = fig.add_subplot(grid[71:100, 0:6])
        plt.plot(ts, s1_t1_spikes, color = s1_t1_c, linewidth = 2, label = 'State A, Cue 1')
        plt.plot(ts, s1_t2_spikes, color = s1_t2_c, linewidth = 2, label = 'State A, Cue 2')
        plt.plot(ts, s2_t1_spikes, color = s2_t1_c, linewidth = 2, label = 'State B, Cue 1')
        plt.plot(ts, s2_t2_spikes, color = s2_t2_c, linewidth = 2, label = 'State B, Cue 2')
        plt.legend(ncol=2, loc = 'upper left', prop={'size': 6})
        #plt.plot([0,0], [2,20], color = 'black', linewidth = 1)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.xlim([-300, 500])
        #plt.ylim([2,20])
        plt.xlabel('Time from State Cue (ms)')
        plt.ylabel('Firing Rate (Hz)')
        fig.suptitle(u_details)
        # save figure
        plt.savefig(savedir + u_details + '.jpg', transparent = True)
        plt.show()

def plot_cue_pop_summary(cue_ttests):

    # which cells are OFC?
    ofc_ix = cue_ttests['ofc'] ==1

    # get proportions of sig units in each area
    n_state_hpc = np.sum(cue_ttests.loc[~ofc_ix, 'state_p'] < .01)
    n_state_ofc = np.sum(cue_ttests.loc[ofc_ix, 'state_p'] < .01)
    n_mix_hpc = np.sum(cue_ttests.loc[~ofc_ix, 'mix_p'] < .01)
    n_mix_ofc = np.sum(cue_ttests.loc[ofc_ix, 'mix_p'] < .01)
    n_hpc = np.sum(~ofc_ix)
    n_ofc = np.sum(ofc_ix)

    # make a contingency table
    hpc_props = np.array([n_state_hpc/n_hpc , n_mix_hpc/n_hpc]) *100
    ofc_props = np.array([n_state_ofc/n_ofc , n_mix_ofc/n_ofc]) *100

    width = 0.35  # the width of the bars
    x = np.arange(2)
    prop_fig = plt.figure(figsize=(6, 6), dpi=300)
    rects1 = plt.bar(x - width/2, hpc_props, width, label='HPC')
    rects2 = plt.bar(x + width/2, ofc_props, width, label='OFC')
    plt.gca().set_xticks([0,1])
    plt.gca().set_xticklabels(['State','Mix'], fontsize =14)
    plt.ylabel('% units with p < .01', fontsize =14)
    plt.xlabel('Factor', fontsize =14)
    plt.legend(['HPC','OFC'], fontsize =14)
    plt.ylim([0, 19])

    stats.chi2_contingency([[n_state_hpc,n_mix_hpc],[n_hpc,n_hpc]],correction=False)

    # define a colormap
    cmap = plt.cm.Paired(np.linspace(0,1,12))


    #-----------------------------------
    # show cohen's D
    hpc_state_mean = np.mean((cue_ttests.loc[~ofc_ix,'state_d']))
    hpc_state_sem = np.std((cue_ttests.loc[~ofc_ix,'state_d'])) / np.sqrt(np.sum(~ofc_ix))

    hpc_mix_mean = np.mean(cue_ttests.loc[~ofc_ix,'mix_d'])
    hpc_mix_sem = np.std(cue_ttests.loc[~ofc_ix,'mix_d']) / np.sqrt(np.sum(~ofc_ix))

    ofc_state_mean = np.mean(np.abs(cue_ttests.loc[ofc_ix,'state_d']))
    ofc_state_sem = np.std(np.abs(cue_ttests.loc[ofc_ix,'state_d'])) / np.sqrt(np.sum(ofc_ix))

    ofc_mix_mean = np.mean(np.abs(cue_ttests.loc[ofc_ix,'mix_d']))
    ofc_mix_sem = np.std(np.abs(cue_ttests.loc[ofc_ix,'mix_d'])) / np.sqrt(np.sum(ofc_ix))


    d_fig = plt.figure(figsize=(6, 6), dpi=300)
    prop_fig = plt.figure(figsize=(6, 6), dpi=300)
    rects1 = plt.bar(x - width/2, [hpc_state_mean,hpc_mix_mean], width,
                    yerr = [hpc_state_sem, hpc_mix_sem],error_kw=dict(lw=2, capsize=0, capthick=0), label='HPC')
    rects2 = plt.bar(x + width/2, [ofc_state_mean,ofc_mix_mean], width, 
                    yerr = [ofc_state_sem, ofc_mix_sem],error_kw=dict(lw=2, capsize=0, capthick=0), label='OFC')
    plt.gca().set_xticks([0,1])
    plt.gca().set_xticklabels(['State','Mix'], fontsize =14)
    plt.ylabel('abs(Cohen\'s D)', fontsize =14)
    plt.xlabel('Factor', fontsize =14)
    plt.legend(['HPC','OFC'], fontsize =14)

    # let's do an ANOVA
    df = pd.DataFrame()
    df['dval'] = np.concatenate((cue_ttests.loc[~ofc_ix,'state_d'], cue_ttests.loc[~ofc_ix,'mix_d'],
                                 cue_ttests.loc[ofc_ix,'state_d'], cue_ttests.loc[ofc_ix,'mix_d']))

    df['brain_area'] = np.concatenate((np.ones(shape = (np.sum(~ofc_ix), )), np.ones(shape = (np.sum(~ofc_ix), )),
                    np.zeros(shape = (np.sum(ofc_ix), )), np.zeros(shape = (np.sum(ofc_ix), ))), axis=0).astype(int)

    df['grouping'] = np.concatenate((np.ones(shape = (np.sum(~ofc_ix), )), np.zeros(shape = (np.sum(~ofc_ix), )),
                    np.ones(shape = (np.sum(ofc_ix), )), np.zeros(shape = (np.sum(ofc_ix), ))), axis=0).astype(int)

    df['unit_id'] = np.concatenate((cue_ttests.loc[~ofc_ix,'cell_id'], cue_ttests.loc[~ofc_ix,'cell_id'],
                                 cue_ttests.loc[ofc_ix,'cell_id'], cue_ttests.loc[ofc_ix,'cell_id']))

    aov = pg.anova(dv = 'dval', between=['brain_area','grouping'], data =df)


def make_cue_psuedo_pop(datadir, n_boots):

    #initialize outputs
    all_train=[None] * n_boots
    all_test=[None] * n_boots
    all_ofc = np.array([])

    # find file names
    datafiles = glob.glob(datadir + '*.mat')

    for fnum, fname in enumerate(datafiles):
            print('file ' + str(fnum+1) + ' \ ' + str(np.size(datafiles)))

            # load this file
            f_data = mat73.loadmat(fname)

            # pull out this file's data
            bhv= pd.DataFrame(f_data['bhv'], columns=f_data['bhv_varnames'])
            cue_spikes = f_data['cue_zFRs']
            ts = f_data['cue_down_times']
            OFC_ix = (f_data['OFC_ix'] == 1).astype(int)
            unit_names = f_data['unit_names']
            # get the mean firing rates for period right after cue-on
            win_start = np.argmin(np.abs(ts - 50))
            win_end = np.argmin(np.abs(ts - 250))

            train_frs = np.mean(cue_spikes[:,win_start:win_end,:], axis=1)

            # delete the big dict to save some space
            del f_data

            file_train = []
            file_test = []
            # collect n_boots random partitions
            for boot in range(n_boots):
                # get the indices of the relevant conditions
                s1_t1_ix =  np.random.permutation(np.where((bhv['state'] == 1) & (bhv['state_type'] == 1))[0])
                s1_t2_ix =  np.random.permutation(np.where((bhv['state'] == 1) & (bhv['state_type'] == -1))[0])
                s2_t1_ix =  np.random.permutation(np.where((bhv['state'] == -1) & (bhv['state_type'] == 1))[0])
                s2_t2_ix =  np.random.permutation(np.where((bhv['state'] == -1) & (bhv['state_type'] == -1))[0])

                # take the first 180 of these random trials for training
                file_train.insert(boot,np.concatenate((train_frs[s1_t1_ix[0:180]],
                                                  train_frs[s1_t2_ix[0:180]],
                                                  train_frs[s2_t1_ix[0:180]],
                                                  train_frs[s1_t2_ix[0:180]]),axis = 0))

                # take the subsequent 10 of these random trials for testing
                file_test.insert(boot,np.concatenate((cue_spikes[s1_t1_ix[190:200],:,:],
                                                 cue_spikes[s1_t2_ix[190:200],:,:],
                                                 cue_spikes[s2_t1_ix[190:200],:,:],
                                                 cue_spikes[s1_t2_ix[190:200],:,:]),axis = 0))

            # now join the splits across files
            if fnum ==0:
                all_train = file_train
                all_test  = file_test
            else:
                for b in range(n_boots):
                    all_train[b] = np.concatenate((all_train[b], file_train[b]), axis=1)
                    all_test[b] = np.concatenate((all_test[b], file_test[b]), axis=2)

            all_ofc = np.concatenate((all_ofc,OFC_ix),axis=0)
    
    train_state_labels = np.concatenate((np.zeros(len(s1_t1_ix[0:180])),
                                         np.zeros(len(s1_t2_ix[0:180])),
                                         np.ones(len(s2_t1_ix[0:180])),
                                         np.ones(len(s2_t2_ix[0:180]))), axis = 0).astype(int)

    test_state_labels = np.concatenate((np.zeros(len(s1_t1_ix[190:200])),
                                        np.zeros(len(s1_t2_ix[190:200])),
                                        np.ones(len(s2_t1_ix[190:200])),
                                        np.ones(len(s2_t2_ix[190:200]))), axis = 0).astype(int)


    return all_train, all_test, train_state_labels, test_state_labels, ts, all_ofc

def cue_within_region_phaselocking(cue_ttests):
    ofc_ix = cue_ttests['ofc'] ==1
    state_ix = cue_ttests['state_p'] < .01
    phase_locked = cue_ttests['rayleigh_p'] < .01

    p_ofc_state_phase_locked = np.sum(phase_locked[ofc_ix & state_ix]) / np.sum(ofc_ix & state_ix)
    p_hpc_state_phase_locked = np.sum(phase_locked[~ofc_ix & state_ix]) / np.sum(~ofc_ix & state_ix)

    p_ofc_Xstate_phase_locked = np.sum(phase_locked[ofc_ix & ~state_ix]) / np.sum(ofc_ix & ~state_ix)
    p_hpc_Xstate_phase_locked = np.sum(phase_locked[~ofc_ix & ~state_ix]) / np.sum(~ofc_ix & ~state_ix)

    # now look at theta modulation across all units
    p_ofc_phase_locked = np.sum(phase_locked[ofc_ix]) / np.sum(ofc_ix)
    p_hpc_phase_locked = np.sum(phase_locked[~ofc_ix]) / np.sum(~ofc_ix)

    hpc_props = np.array([p_hpc_state_phase_locked , p_hpc_Xstate_phase_locked, p_hpc_phase_locked]) *100
    ofc_props = np.array([p_ofc_state_phase_locked , p_ofc_Xstate_phase_locked, p_ofc_phase_locked]) *100

    width = 0.35  # the width of the bars
    x = np.arange(3)
    prop_fig = plt.figure(figsize=(6, 6), dpi=300)
    rects1 = plt.bar(x - width/2, hpc_props, width, label='HPC')
    rects2 = plt.bar(x + width/2, ofc_props, width, label='OFC')
    plt.gca().set_xticks([0,1,2])
    plt.gca().set_xticklabels(['state p<.01','state p>.01','all units'], fontsize =14)
    plt.ylabel('% phase-modulated units (p<.01)', fontsize =14)
    plt.xlabel('Grouping', fontsize =14)
    plt.legend(['HPC','OFC'], fontsize =14)
    plt.title('${\Theta}$ modulates state-coding units during cue phase')
    plt.savefig('cue_phase_locking.svg', transparent = True)


def choice_within_region_phaselocking(choice_glms):
    ofc_ix = choice_glms['ofc'] ==1
    state_ix = (choice_glms['state_p'] < .01) | (choice_glms['state:val_p'] < .01)
    phase_locked = choice_glms['rayleigh_p'] < .01

    # look at theta modulation across all units
    p_hpc_phase_locked = np.sum(phase_locked[~ofc_ix]) / np.sum(~ofc_ix)
    p_ofc_phase_locked = np.sum(phase_locked[ofc_ix]) / np.sum(ofc_ix)

    p_hpc_state_phase_locked = np.sum(phase_locked[~ofc_ix & state_ix]) / np.sum(~ofc_ix & state_ix)
    p_ofc_state_phase_locked = np.sum(phase_locked[ofc_ix & state_ix]) / np.sum(ofc_ix & state_ix)

    p_hpc_Xstate_phase_locked = np.sum(phase_locked[~ofc_ix & ~state_ix]) / np.sum(~ofc_ix & ~state_ix)
    p_ofc_Xstate_phase_locked = np.sum(phase_locked[ofc_ix & ~state_ix]) / np.sum(ofc_ix & ~state_ix)



    hpc_props = np.array([p_hpc_state_phase_locked , p_hpc_Xstate_phase_locked, p_hpc_phase_locked]) *100
    ofc_props = np.array([p_ofc_state_phase_locked , p_ofc_Xstate_phase_locked, p_ofc_phase_locked]) *100

    width = 0.35  # the width of the bars
    x = np.arange(3)
    prop_fig = plt.figure(figsize=(6, 6), dpi=300)
    rects1 = plt.bar(x - width/2, hpc_props, width, label='HPC')
    rects2 = plt.bar(x + width/2, ofc_props, width, label='OFC')
    plt.gca().set_xticks([0,1,2])
    plt.gca().set_xticklabels(['state p<.01','state p>.01','all units'], fontsize =14)
    plt.ylabel('% phase-modulated units (p<.01)', fontsize =14)
    plt.xlabel('Grouping', fontsize =14)
    plt.legend(['HPC','OFC'], fontsize =14)
    plt.title('${\Theta}$ modulates state-coding units during choice phase')
    plt.savefig('choice_phase_locking.svg', transparent = True)


def decode_cue_psuedo_pop(all_train, all_test, train_state_labels, test_state_labels, ts, all_ofc):

    n_boots = len(all_train)
    n_times = len(ts)

    # initialize output
    hpc_state_acc = np.zeros(shape = (n_boots, n_times))
    ofc_state_acc = np.zeros(shape = (n_boots, n_times))
    hpc_mix_acc = np.zeros(shape = (n_boots, n_times))
    ofc_mix_acc = np.zeros(shape = (n_boots, n_times))

    hpc_mdl = LinearDiscriminantAnalysis()
    ofc_mdl = LinearDiscriminantAnalysis()


    for b in range(n_boots):

        # pull out this boot's firing rates
        hpc_mdl.fit(all_train[b][:,all_ofc==0], train_state_labels)
        ofc_mdl.fit(all_train[b][:,all_ofc==1], train_state_labels)

        for t in range(n_times):
            print('fold ' + str(b) + ', ts ' + str(t))

            hpc_state_acc[b,t] = np.mean(hpc_mdl.predict(all_test[b][:,t,all_ofc==0]) == test_state_labels)
            ofc_state_acc[b,t] = np.mean(ofc_mdl.predict(all_test[b][:,t,all_ofc==1]) == test_state_labels)


    cmap = plt.cm.Paired(np.linspace(0,1,12))
    plt.plot(ts,np.mean(hpc_state_acc,axis=0),color=cmap[5,:], label = 'HPC state')
    plt.plot(ts,np.mean(ofc_state_acc,axis=0),color=cmap[1,:], label = 'OFC state')
    plt.legend()
    plt.ylabel('accuracy')
    plt.xlabel('time from state cue on (ms)')

def plot_choice_selective_units(choice_glms):

    savedir = 'C:/Users/Thomas Elston/Documents/PYTHON/CuedFlexVal/choice_units/'
    # plot rasters and psth aligned to cue on
    sig_coding_ix = (choice_glms['state_p'] < 1) | (choice_glms['val_p'] < 1)
    sig_cells = choice_glms.loc[sig_coding_ix,:]
    sig_cells = sig_cells.reset_index(drop=True)

    # loop over individual cells and plot
    for u in range(len(sig_cells)):

        print(str(u) +  ' / ' + str(len(sig_cells)))

        # load this unit's file and extract this unit
        f_data = mat73.loadmat(sig_cells['file'][u])
        spikes = f_data['raw_choice_FRs']
        rasters = f_data['choice_rasters']
        u_name = sig_cells['u_name'][u][0]
        u_in_ofc = sig_cells['ofc'][u]==1
        ts = f_data['choice_raw_times']
        ts = ts[0:len(rasters[0,:,0])]
        bhv= pd.DataFrame(f_data['bhv'], columns=f_data['bhv_varnames'])
        del f_data

        trials2use = (bhv['pickedbest']==1) 
        bhv = bhv.loc[trials2use,:]
        spikes = spikes[trials2use,:,:]
        rasters = rasters[trials2use,:,:]

        # pull out this unit's spikes and rasters
        u_spikes = spikes[:,:,(sig_cells['unit_num'][u]).astype(int)]
        u_rasters = rasters[:,:,(sig_cells['unit_num'][u]).astype(int)]

        # do a sliding window glm for the period around the choice
        sig_factors = sliding_window_glm(bhv, u_spikes, ts, 50)

        # get indices for the states and cues
        s1_ix = bhv['state'] ==1
        s2_ix = bhv['state'] ==-1
        t1_ix = bhv['state_type'] ==1
        t2_ix = bhv['state_type'] ==-1

        # get indices for values
        v1_ix = bhv['chosenval'] ==1
        v2_ix = bhv['chosenval'] ==2
        v3_ix = bhv['chosenval'] ==3
        v4_ix = bhv['chosenval'] ==4

        # get state-type combination spikes
        s1_t1_spikes = np.mean(u_spikes[s1_ix & t1_ix,:],axis =0)
        s1_t2_spikes = np.mean(u_spikes[s1_ix & t2_ix,:],axis =0)
        s2_t1_spikes = np.mean(u_spikes[s2_ix & t1_ix,:],axis =0)
        s2_t2_spikes = np.mean(u_spikes[s2_ix & t2_ix,:],axis =0)

        # get value spikes
        v1_spikes = np.mean(u_spikes[v1_ix,:],axis =0)
        v2_spikes = np.mean(u_spikes[v2_ix,:],axis =0)
        v3_spikes = np.mean(u_spikes[v3_ix,:],axis =0)
        v4_spikes = np.mean(u_spikes[v4_ix,:],axis =0)

        # n_trials to show rasters for each condition
        n_trials = 200
        n_val_trials = 60;
        # grab a random n_trials for raster plotting
        s1_t1_rasters = u_rasters[np.random.permutation(np.where(s1_ix & t1_ix)[0])[0:n_trials],:]
        s1_t2_rasters = u_rasters[np.random.permutation(np.where(s1_ix & t2_ix)[0])[0:n_trials],:]
        s2_t1_rasters = u_rasters[np.random.permutation(np.where(s2_ix & t1_ix)[0])[0:n_trials],:]
        s2_t2_rasters = u_rasters[np.random.permutation(np.where(s2_ix & t2_ix)[0])[0:n_trials],:]

        v1_rasters = u_rasters[np.random.permutation(np.where(v1_ix)[0])[0:n_val_trials],:]
        v2_rasters = u_rasters[np.random.permutation(np.where(v2_ix)[0])[0:n_val_trials],:]
        v3_rasters = u_rasters[np.random.permutation(np.where(v3_ix)[0])[0:n_val_trials],:]
        v4_rasters = u_rasters[np.random.permutation(np.where(v4_ix)[0])[0:n_val_trials],:]

        # define colormaps
        state_cmap = plt.cm.Paired(np.linspace(0,1,12))
        # define the colors for the states/types
        s1_t1_c = state_cmap[1,:] # sat blue
        s1_t2_c = state_cmap[0,:] # unsat blue
        s2_t1_c = state_cmap[5,:] # sat red
        s2_t2_c = state_cmap[4,:] # unsat red

        val_cmap = plt.cm.tab20b(np.linspace(0,1,20))
        v1_c = val_cmap[3,:]
        v2_c = val_cmap[2,:]
        v3_c = val_cmap[1,:]
        v4_c = val_cmap[0,:]

        # put the rasters into a big array and make a corresponding array of colors
        state_rasters = np.concatenate((s1_t1_rasters,s1_t2_rasters,s2_t1_rasters,s2_t2_rasters), axis=0)
        state_colors  = np.concatenate((np.matlib.repmat(s1_t1_c,n_trials,1),
                                      np.matlib.repmat(s1_t2_c,n_trials,1),
                                      np.matlib.repmat(s2_t1_c,n_trials,1),
                                      np.matlib.repmat(s2_t2_c,n_trials,1)), axis = 0)

        val_rasters = np.concatenate((v1_rasters,v2_rasters,v3_rasters,v4_rasters), axis=0)
        val_colors  = np.concatenate((np.matlib.repmat(v1_c,n_val_trials,1),
                                      np.matlib.repmat(v2_c,n_val_trials,1),
                                      np.matlib.repmat(v3_c,n_val_trials,1),
                                      np.matlib.repmat(v4_c,n_val_trials,1)), axis = 0)                              

        # update the row values of the rasters so we can plot row by row  
        for t in range(len(state_rasters)): 
            state_rasters[t,:] = (state_rasters[t,:]*t) 
        state_rasters[state_rasters == 0] = np.nan   

        for t in range(len(val_rasters)): 
            val_rasters[t,:] = (val_rasters[t,:]*t) 
        val_rasters[val_rasters == 0] = np.nan   

        # get details of unit name
        if u_in_ofc:
            u_details = 'OFC_' + u_name
        else:
            u_details = 'HPC_' + u_name


        #---------------------------------------
        #               PLOT
        #---------------------------------------
        sns.set_style("ticks")
        fig = plt.figure(figsize=(8, 8), dpi=300)
        grid = fig.add_gridspec(100, 21)
        # plot the rasters
        raster_ax = fig.add_subplot(grid[0:60, 0:9])
        for t in range(len(state_rasters)): 
            plt.scatter(ts,state_rasters[t,:], color = state_colors[t,:], s=.2)
        plt.gca().axis('off')
        plt.xlim([-1000, 1000])
        lw = 2

        psth_ax = fig.add_subplot(grid[61:100, 0:9])
        plt.plot(ts, s1_t1_spikes, color = s1_t1_c, linewidth = lw, label = 'State A, Cue 1')
        plt.plot(ts, s1_t2_spikes, color = s1_t2_c, linewidth = lw, label = 'State A, Cue 2')
        plt.plot(ts, s2_t1_spikes, color = s2_t1_c, linewidth = lw, label = 'State B, Cue 1')
        plt.plot(ts, s2_t2_spikes, color = s2_t2_c, linewidth = lw, label = 'State B, Cue 2')
        plt.plot(ts, sig_factors[:,0]*12,color = 'k', linewidth = 4)
        #plt.legend(ncol=2, loc = 'upper left', prop={'size': 6})
        #plt.plot([-700,-700], [1,8], color = 'black', linewidth = 1)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.xlim([-1000, 1000])
        
        #plt.ylim([10,50])
        plt.xlabel('Time from Pics On (ms)')
        plt.ylabel('Firing Rate (Hz)')

        # VALUE
        # plot the value rasters
        val_raster_ax = fig.add_subplot(grid[0:60, 11:21])
        for t in range(len(val_rasters)): 
            plt.scatter(ts,val_rasters[t,:], color = val_colors[t,:], s=.2)
        plt.gca().axis('off')
        plt.xlim([-1000, 1000])
        #plt.plot([0,0], plt.gca().get_ylim(), color = 'black')

        psth_ax = fig.add_subplot(grid[61:100, 11:21])
        plt.plot(ts, v1_spikes, color = v1_c, linewidth = lw, label = 'val = 1')
        plt.plot(ts, v2_spikes, color = v2_c, linewidth = lw, label = 'val = 2')
        plt.plot(ts, v3_spikes, color = v3_c, linewidth = lw, label = 'val = 3')
        plt.plot(ts, v4_spikes, color = v4_c, linewidth = lw, label = 'val = 4')
        plt.plot(ts, sig_factors[:,1]*12,color = 'k', linewidth = 4)
       # plt.legend(ncol=2, loc = 'upper left', prop={'size': 6})
        #plt.plot([0,0], [2,20], color = 'black', linewidth = 1)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.xlim([-1000, 1000])
        #plt.ylim([10,50])
        
        #plt.xlabel('Time from Pics On (ms)')
        #plt.ylabel('Firing Rate (Hz)')
        fig.suptitle(u_details)
        # save figure
        #plt.savefig(savedir + u_details + '.jpg', transparent = True)
        #plt.savefig(u_details + '.svg', transparent = True)

        plt.show()


def sliding_window_glm(bhv, u_spikes, ts, time_thresh):

    n_trials, n_times = np.shape(u_spikes)
    sig_factors = np.zeros(shape=(n_times, 2))

    # only assess times around the choice
    win_start = np.argmin(np.abs(ts - -700))
    win_end = np.argmin(np.abs(ts - 1000)) 
    win_times = np.arange(win_start,win_end)
    glm_df = pd.DataFrame()

    glm_df['state'] = bhv['state']
    glm_df['val'] = bhv['chosenval']
    glm_df['t_num'] = np.arange(n_trials)

    for w in win_times:
        glm_df['fr'] = u_spikes[:,w]

        mdl = smf.glm(formula = 'fr ~ state*val + t_num',data = glm_df).fit()
        sig_factors[w,0] = mdl.pvalues['state'] < .05
        sig_factors[w,1] = mdl.pvalues['val'] < .05

    state_lens, state_pos, state_seq_type = ut.find_sequences(sig_factors[:,0])
    val_lens, val_pos, val_seq_type = ut.find_sequences(sig_factors[:,1])

    state_pos = state_pos[(state_seq_type ==1) &  (state_lens > time_thresh)]
    state_lens = state_lens[(state_seq_type ==1) & (state_lens > time_thresh)]
    val_pos = val_pos[(val_seq_type ==1) &  (val_lens > time_thresh)]
    val_lens = val_lens[(val_seq_type ==1) & (val_lens > time_thresh)]

    new_sig_factors = np.zeros_like(sig_factors)

    # now fill out new_sig_factors
    for sr in range(len(state_pos)):
        new_sig_factors[state_pos[sr]:state_pos[sr]+state_lens[sr],0]=1

    for vr in range(len(val_pos)):
        new_sig_factors[val_pos[vr]:val_pos[vr]+val_lens[vr],1]=1

    new_sig_factors[new_sig_factors==0] = np.nan

    return new_sig_factors

def plot_choice_pop_summary(choice_glms):

    # which cells are OFC?
    ofc_ix = choice_glms['ofc'] ==1

    # get proportions of sig units in each area
    thresh = .01
    state_ix = (choice_glms['state_p'] < thresh) & (choice_glms['state:val_p'] > thresh)
    val_ix = (choice_glms['val_p'] < thresh)  & (choice_glms['state:val_p'] > thresh)
    state_x_val_ix = (choice_glms['state:val_p'] < thresh)

    n_state_hpc = np.sum(state_ix & ~ofc_ix)
    n_state_ofc = np.sum(state_ix & ofc_ix)
    n_val_hpc = np.sum(val_ix & ~ofc_ix)
    n_val_ofc = np.sum(val_ix & ofc_ix)
    n_state_val_hpc = np.sum(state_x_val_ix & ~ofc_ix)
    n_state_val_ofc = np.sum(state_x_val_ix & ofc_ix)
    n_hpc = np.sum(~ofc_ix)
    n_ofc = np.sum(ofc_ix)

    # make a contingency table
    hpc_props = np.array([n_state_hpc/n_hpc , n_val_hpc/n_hpc, n_state_val_hpc/n_hpc]) *100
    ofc_props = np.array([n_state_ofc/n_ofc , n_val_ofc/n_ofc, n_state_val_ofc/n_ofc]) *100

    width = 0.35  # the width of the bars
    x = np.arange(3)
    prop_fig = plt.figure(figsize=(6, 6), dpi=300)
    rects1 = plt.bar(x - width/2, hpc_props, width, label='HPC')
    rects2 = plt.bar(x + width/2, ofc_props, width, label='OFC')
    plt.gca().set_xticks([0,1,2])
    plt.gca().set_xticklabels(['State only','Value only','State x Value'], fontsize =14)
    plt.ylabel('% units with p < .01', fontsize =14)
    plt.xlabel('Factor', fontsize =14)
    plt.legend(['HPC','OFC'], fontsize =14)
    #plt.ylim([0, 25])
    plt.savefig('choice_props.svg', transparent = True)


    #stats.chi2_contingency([[n_val_hpc,n_val_ofc],[n_hpc,n_ofc]],correction=False)



    #-----------------------------------
    # show beta weights
    hpc_state_mean = np.mean(np.abs(choice_glms.loc[~ofc_ix,'state_b']))
    hpc_state_sem = np.std(np.abs(choice_glms.loc[~ofc_ix,'state_b'])) / np.sqrt(np.sum(~ofc_ix))

    hpc_val_mean = np.mean(np.abs(choice_glms.loc[~ofc_ix,'val_b']))
    hpc_val_sem = np.std(np.abs(choice_glms.loc[~ofc_ix,'val_b'])) / np.sqrt(np.sum(~ofc_ix))

    ofc_state_mean = np.mean(np.abs(choice_glms.loc[ofc_ix,'state_b']))
    ofc_state_sem = np.std(np.abs(choice_glms.loc[ofc_ix,'state_b'])) / np.sqrt(np.sum(ofc_ix))

    ofc_val_mean = np.mean(np.abs(choice_glms.loc[ofc_ix,'val_b']))
    ofc_val_sem = np.std(np.abs(choice_glms.loc[ofc_ix,'val_b'])) / np.sqrt(np.sum(ofc_ix))


    b_fig = plt.figure(figsize=(6, 6), dpi=300)
    prop_fig = plt.figure(figsize=(6, 6), dpi=300)
    rects1 = plt.bar(x - width/2, [hpc_state_mean,hpc_val_mean], width,
                    yerr = [hpc_state_sem, hpc_val_sem],error_kw=dict(lw=2, capsize=0, capthick=0), label='HPC')
    rects2 = plt.bar(x + width/2, [ofc_state_mean,ofc_val_mean], width, 
                    yerr = [ofc_state_sem, ofc_val_sem],error_kw=dict(lw=2, capsize=0, capthick=0), label='OFC')
    plt.gca().set_xticks([0,1])
    plt.gca().set_xticklabels(['State','Value'], fontsize =14)
    plt.ylabel('abs(\u03B2)', fontsize =14)
    plt.xlabel('Factor', fontsize =14)
    plt.legend(['HPC','OFC'], fontsize =14)

def make_choice_psuedo_pop(datadir, n_boots):

    #initialize outputs
    all_train=[None] * n_boots
    all_test=[None] * n_boots
    all_ofc = np.array([])

    # find file names
    datafiles = glob.glob(datadir + '*.mat')

    for fnum, fname in enumerate(datafiles):
            print('file ' + str(fnum+1) + ' \ ' + str(np.size(datafiles)))

            # load this file
            f_data = mat73.loadmat(fname)

            # pull out this file's data
            bhv= pd.DataFrame(f_data['bhv'], columns=f_data['bhv_varnames'])
            choice_spikes = f_data['down_choice_FRs']
            ts = f_data['choice_down_times']
            OFC_ix = (f_data['OFC_ix'] == 1).astype(int)
            unit_names = f_data['unit_names']
            # get the mean firing rates for period right after cue-on
            win_start = np.argmin(np.abs(ts - 50))
            win_end = np.argmin(np.abs(ts - 250))

            train_frs = np.mean(choice_spikes[:,win_start:win_end,:], axis=1)

            # delete the big dict to save some space
            del f_data

            file_train = []
            file_test = []
            # collect n_boots random partitions
            for boot in range(n_boots):
                # get the indices of the relevant conditions
                s1_t1_ix =  np.random.permutation(np.where((bhv['state'] == 1) & (bhv['state_type'] == 1))[0])
                s1_t2_ix =  np.random.permutation(np.where((bhv['state'] == 1) & (bhv['state_type'] == -1))[0])
                s2_t1_ix =  np.random.permutation(np.where((bhv['state'] == -1) & (bhv['state_type'] == 1))[0])
                s2_t2_ix =  np.random.permutation(np.where((bhv['state'] == -1) & (bhv['state_type'] == -1))[0])

                # take the first 180 of these random trials for training
                file_train.insert(boot,np.concatenate((train_frs[s1_t1_ix[0:180]],
                                                  train_frs[s1_t2_ix[0:180]],
                                                  train_frs[s2_t1_ix[0:180]],
                                                  train_frs[s1_t2_ix[0:180]]),axis = 0))

                # take the subsequent 10 of these random trials for testing
                file_test.insert(boot,np.concatenate((choice_spikes[s1_t1_ix[190:200],:,:],
                                                 choice_spikes[s1_t2_ix[190:200],:,:],
                                                 choice_spikes[s2_t1_ix[190:200],:,:],
                                                 choice_spikes[s1_t2_ix[190:200],:,:]),axis = 0))

            # now join the splits across files
            if fnum ==0:
                all_train = file_train
                all_test  = file_test
            else:
                for b in range(n_boots):
                    all_train[b] = np.concatenate((all_train[b], file_train[b]), axis=1)
                    all_test[b] = np.concatenate((all_test[b], file_test[b]), axis=2)

            all_ofc = np.concatenate((all_ofc,OFC_ix),axis=0)
    
    train_state_labels = np.concatenate((np.zeros(len(s1_t1_ix[0:180])),
                                         np.zeros(len(s1_t2_ix[0:180])),
                                         np.ones(len(s2_t1_ix[0:180])),
                                         np.ones(len(s2_t2_ix[0:180]))), axis = 0).astype(int)

    test_state_labels = np.concatenate((np.zeros(len(s1_t1_ix[190:200])),
                                        np.zeros(len(s1_t2_ix[190:200])),
                                        np.ones(len(s2_t1_ix[190:200])),
                                        np.ones(len(s2_t2_ix[190:200]))), axis = 0).astype(int)


    return all_train, all_test, train_state_labels, test_state_labels, ts, all_ofc


def make_choice_value_psuedo_pop(datadir, n_boots):

    #initialize outputs
    all_train=[None] * n_boots
    all_test=[None] * n_boots
    all_ofc = np.array([])

    # find file names
    datafiles = glob.glob(datadir + '*.mat')

    for fnum, fname in enumerate(datafiles):
            print('file ' + str(fnum+1) + ' \ ' + str(np.size(datafiles)))

            # load this file
            f_data = mat73.loadmat(fname)

            # pull out this file's data
            bhv= pd.DataFrame(f_data['bhv'], columns=f_data['bhv_varnames'])
            choice_spikes = f_data['down_choice_FRs']
            ts = f_data['choice_down_times']
            OFC_ix = (f_data['OFC_ix'] == 1).astype(int)
            unit_names = f_data['unit_names']
            # get the mean firing rates for period right after cue-on
            win_start = np.argmin(np.abs(ts - 50))
            win_end = np.argmin(np.abs(ts - 250))

            train_frs = np.mean(choice_spikes[:,win_start:win_end,:], axis=1)

            # delete the big dict to save some space
            del f_data

            file_train = []
            file_test = []
            # collect n_boots random partitions
            for boot in range(n_boots):
                # get the indices of the relevant conditions
                val1_ix =  np.random.permutation(np.where((bhv['chosenval'] == 1))[0])
                val2_ix =  np.random.permutation(np.where((bhv['chosenval'] == 2))[0])
                val3_ix =  np.random.permutation(np.where((bhv['chosenval'] == 3))[0])
                val4_ix =  np.random.permutation(np.where((bhv['chosenval'] == 4))[0])


                # # take the first 180 of these random trials for training
                # file_train.insert(boot,np.concatenate((train_frs[val1_ix[0:50]],
                #                                   train_frs[val2_ix[0:50]],
                #                                   train_frs[val3_ix[0:50]],
                #                                   train_frs[val4_ix[0:50]]),axis = 0))

                # # take the subsequent 10 of these random trials for testing
                # file_test.insert(boot,np.concatenate((choice_spikes[val1_ix[-1 - 10:-1],:,:],
                #                                  choice_spikes[val2_ix[-1 - 10:-1],:,:],
                #                                  choice_spikes[val3_ix[-1 - 10:-1],:,:],
                #                                  choice_spikes[val4_ix[-1 - 10:-1],:,:]),axis = 0))

                                # take the first 180 of these random trials for training
                file_train.insert(boot,np.concatenate((train_frs[val2_ix[0:50]],
                                                       train_frs[val3_ix[0:50]],
                                                       train_frs[val4_ix[0:50]]),axis = 0))

                # take the subsequent 10 of these random trials for testing
                file_test.insert(boot,np.concatenate((choice_spikes[val2_ix[-1 - 10:-1],:,:],
                                                 choice_spikes[val3_ix[-1 - 10:-1],:,:],
                                                 choice_spikes[val4_ix[-1 - 10:-1],:,:]),axis = 0))

            # now join the splits across files
            if fnum ==0:
                all_train = file_train
                all_test  = file_test
            else:
                for b in range(n_boots):
                    all_train[b] = np.concatenate((all_train[b], file_train[b]), axis=1)
                    all_test[b] = np.concatenate((all_test[b], file_test[b]), axis=2)

            all_ofc = np.concatenate((all_ofc,OFC_ix),axis=0)
    
    # train_val_labels = np.concatenate((np.ones(len(val1_ix[0:50])),
    #                                     np.ones(len(val2_ix[0:50]))*2,
    #                                     np.ones(len(val3_ix[0:50]))*3,
    #                                     np.ones(len(val4_ix[0:50]))*4), axis=0).astype(int)

    # test_val_labels = np.concatenate((np.ones(len(val1_ix[-1 - 10:-1])),
    #                                     np.ones(len(val2_ix[-1 - 10:-1]))*2,
    #                                     np.ones(len(val3_ix[-1 - 10:-1]))*3,
    #                                     np.ones(len(val4_ix[-1 - 10:-1]))*4), axis = 0).astype(int)
    train_val_labels = np.concatenate((np.ones(len(val2_ix[0:50]))*2,
                                        np.ones(len(val3_ix[0:50]))*3,
                                        np.ones(len(val4_ix[0:50]))*4), axis=0).astype(int)

    test_val_labels = np.concatenate((np.ones(len(val2_ix[-1 - 10:-1]))*2,
                                        np.ones(len(val3_ix[-1 - 10:-1]))*3,
                                        np.ones(len(val4_ix[-1 - 10:-1]))*4), axis = 0).astype(int)


    return all_train, all_test, train_val_labels, test_val_labels, ts, all_ofc


def decode_choice_val_psuedo_pop(all_train, all_test, train_val_labels, test_val_labels, ts, all_ofc):

    n_boots = len(all_train)
    n_times = len(ts)

    # initialize output
    hpc_val_acc = np.zeros(shape = (n_boots, n_times))
    ofc_val_acc = np.zeros(shape = (n_boots, n_times))

    hpc_mdl = LinearDiscriminantAnalysis()
    ofc_mdl = LinearDiscriminantAnalysis()


    for b in range(n_boots):

        # pull out this boot's firing rates
        hpc_mdl.fit(all_train[b][:,all_ofc==0], train_val_labels)
        ofc_mdl.fit(all_train[b][:,all_ofc==1], train_val_labels)

        for t in range(n_times):
            print('fold ' + str(b) + ', ts ' + str(t))

            hpc_val_acc[b,t] = np.mean(hpc_mdl.predict(all_test[b][:,t,all_ofc==0]) == test_val_labels)
            ofc_val_acc[b,t] = np.mean(ofc_mdl.predict(all_test[b][:,t,all_ofc==1]) == test_val_labels)


    cmap = plt.cm.Paired(np.linspace(0,1,12))
    plt.plot(ts,np.mean(hpc_val_acc,axis=0),color=cmap[5,:], label = 'HPC val')
    plt.plot(ts,np.mean(ofc_val_acc,axis=0),color=cmap[1,:], label = 'OFC val')
    plt.legend()
    plt.ylabel('accuracy')
    plt.xlabel('time from state choice on (ms)')

def check_val_by_state(choice_glms):

    ofc_ix = choice_glms['ofc']==1
    s1_sig = (choice_glms['s1_val_p']<.01) & (choice_glms['s2_val_p']>.01)
    s2_sig = (choice_glms['s2_val_p']<.01) & (choice_glms['s1_val_p']>.01)
    both_sig = (choice_glms['s1_val_p']<.01) & (choice_glms['s2_val_p']<.01) 

    n_nested = np.sum(s1_sig) + np.sum(s2_sig)
    n_both = np.sum(both_sig)

    sns.set_style("ticks")
    fig = plt.figure(figsize=(8, 8), dpi=300)
    grid = fig.add_gridspec(21, 21)
    # plot the rasters

    hpc_beta_ax = fig.add_subplot(grid[0:9, 0:9])
    plt.scatter(np.abs(choice_glms['s1_val_b'].loc[s1_sig&~ofc_ix]), np.abs(choice_glms['s2_val_b'].loc[s1_sig&~ofc_ix]),label='S1 only')
    plt.scatter(np.abs(choice_glms['s1_val_b'].loc[s2_sig&~ofc_ix]), np.abs(choice_glms['s2_val_b'].loc[s2_sig&~ofc_ix]),label='S2 only')
    plt.scatter(np.abs(choice_glms['s1_val_b'].loc[both_sig&~ofc_ix]), np.abs(choice_glms['s2_val_b'].loc[both_sig&~ofc_ix]),label='both')
    plt.xlim([0,.3])
    plt.ylim([0,.3])
    plt.xlabel('|$\\beta$val| in State 1')
    plt.ylabel('|$\\beta$val| in State 2')
    plt.plot(plt.gca().get_xlim(),plt.gca().get_ylim(),color='k', linestyle = 'dotted', label = 'unity')
    plt.legend()
    plt.title('HPC')
   
    ofc_beta_ax = fig.add_subplot(grid[0:9, 11:21])
    plt.scatter(np.abs(choice_glms['s1_val_b'].loc[s1_sig&ofc_ix]), np.abs(choice_glms['s2_val_b'].loc[s1_sig&ofc_ix]),label='S1 only')
    plt.scatter(np.abs(choice_glms['s1_val_b'].loc[s2_sig&ofc_ix]), np.abs(choice_glms['s2_val_b'].loc[s2_sig&ofc_ix]),label='S2 only')
    plt.scatter(np.abs(choice_glms['s1_val_b'].loc[both_sig&ofc_ix]), np.abs(choice_glms['s2_val_b'].loc[both_sig&ofc_ix]),label='both')
    plt.xlim([0,.3])
    plt.ylim([0,.3])
    plt.plot(plt.gca().get_xlim(),plt.gca().get_ylim(),color='k', linestyle = 'dotted', label = 'unity')
    plt.title('OFC')


    hpc_int_ax = fig.add_subplot(grid[11:21, 0:9])
    plt.scatter(np.abs(choice_glms['s1_val_int'].loc[s1_sig&~ofc_ix]), np.abs(choice_glms['s2_val_int'].loc[s1_sig&~ofc_ix]),label='S1 only')
    plt.scatter(np.abs(choice_glms['s1_val_int'].loc[s2_sig&~ofc_ix]), np.abs(choice_glms['s2_val_int'].loc[s2_sig&~ofc_ix]),label='S2 only')
    plt.scatter(np.abs(choice_glms['s1_val_int'].loc[both_sig&~ofc_ix]), np.abs(choice_glms['s2_val_int'].loc[both_sig&~ofc_ix]),label='both')
    plt.xlim([0,1.4])
    plt.ylim([0,1.4])
    plt.xlabel('INTval in State 1')
    plt.ylabel('INTval in State 2')
    plt.plot(plt.gca().get_xlim(),plt.gca().get_ylim(),color='k', linestyle = 'dotted', label = 'unity')
   
    ofc_int_ax = fig.add_subplot(grid[11:21, 11:21])
    plt.scatter(np.abs(choice_glms['s1_val_int'].loc[s1_sig&ofc_ix]), np.abs(choice_glms['s2_val_int'].loc[s1_sig&ofc_ix]),label='S1 only')
    plt.scatter(np.abs(choice_glms['s1_val_int'].loc[s2_sig&ofc_ix]), np.abs(choice_glms['s2_val_int'].loc[s2_sig&ofc_ix]),label='S2 only')
    plt.scatter(np.abs(choice_glms['s1_val_int'].loc[both_sig&ofc_ix]), np.abs(choice_glms['s2_val_int'].loc[both_sig&ofc_ix]),label='both')
    plt.xlim([0,1.4])
    plt.ylim([0,1.4])
    plt.plot(plt.gca().get_xlim(),plt.gca().get_ylim(),color='k', linestyle = 'dotted', label = 'unity')





















