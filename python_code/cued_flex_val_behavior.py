# -*- coding: utf-8 -*-
"""
cued_flex_val_behavior.py
Loads and parses H5 files generated from NIMH monkeylogic
@author: Thomas Elston
"""

import os
import pandas as pd
import numpy as np
import pingouin as pg
import h5py
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import statsmodels.api as sm
import scipy.signal as sig
import scipy.stats as stats

import utils as ut


def getMLH5(datadir, save_dir = None):
    
    
    # set the path so python is pointed at the data folder
    os.chdir(datadir) 
    
    fnames = os.listdir(datadir)
    
    alldata = pd.DataFrame()
    print('Loading file#: ',end=' ')
    for i in range(len(fnames)):
        
        print(str(i) ,end=' ')
                
        f = h5py.File(fnames[i],'r')
        
        ftrials  = list(f['ML'].keys())[1:-1]        
        sessiondf = pd.DataFrame()
      
        for t in range(len(ftrials)):
                    
            sessiondf.loc[t,1] = fnames[i]
            sessiondf.loc[t,2] = t
            sessiondf.loc[t,3] = f['ML'][ftrials[t]]['UserVars']['UseTrial'][0]
            sessiondf.loc[t,4] = f['ML'][ftrials[t]]['UserVars']['rule'][0]
            sessiondf.loc[t,5] = f['ML'][ftrials[t]]['UserVars']['rtype'][0]
            sessiondf.loc[t,6] = f['ML'][ftrials[t]]['UserVars']['forced'][0]
            
            if sessiondf.loc[t,6] == 0:
            
                sessiondf.loc[t,7] = f['ML'][ftrials[t]]['UserVars']['OptionVals'][0]
                sessiondf.loc[t,8] = f['ML'][ftrials[t]]['UserVars']['OptionVals'][1]
            else:
                if f["ML"][ftrials[t]]['UserVars']['SideChosen'][0] == -1:
                    sessiondf.loc[t,7] = f['ML'][ftrials[t]]['UserVars']['OptionVals'][0]
                    sessiondf.loc[t,8] = np.NaN
                else:
                    sessiondf.loc[t,8] = f['ML'][ftrials[t]]['UserVars']['OptionVals'][0]
                    sessiondf.loc[t,7] = np.NaN
                    
    
            sessiondf.loc[t,9] = f['ML'][ftrials[t]]['UserVars']['ChosenVal'][0]
            sessiondf.loc[t,10] = f['ML'][ftrials[t]]['UserVars']['PickedBestOpt'][0]
            sessiondf.loc[t,11] = f['ML'][ftrials[t]]['UserVars']['RT'][0]
            sessiondf.loc[t,12] = f['ML'][ftrials[t]]['UserVars']['SideChosen'][0]   
            
            
            # get the event code for stim on for saccade detection
            event_codes = np.array([f['ML'][ftrials[t]]['BehavioralCodes']['CodeNumbers']])
            event_times = np.array([f['ML'][ftrials[t]]['BehavioralCodes']['CodeTimes']])
            stim_on_time = np.round(event_times[event_codes==40] / 2).astype(int)
            stim_off_time = np.round(event_times[event_codes==41] / 2).astype(int)
 
            # get eye data for the trial
            eye = np.squeeze(np.array([f['ML'][ftrials[t]]['AnalogData']['Eye']]))
            
            # set some dummy vars for the number + locations of saccades
            sessiondf.loc[t,13] = np.NaN # nsaccs
            sessiondf.loc[t,14] = np.NaN # sacc 1
            sessiondf.loc[t,15] = np.NaN # sacc 2
            sessiondf.loc[t,16] = np.NaN # sacc 3
            sessiondf.loc[t,17] = np.NaN # sacc 4
            sessiondf.loc[t,18] = np.NaN # sacc 5
            sessiondf.loc[t,19] = np.NaN # sacc 1 value

            
            if (stim_on_time.size > 0) & (stim_off_time.size > 0):
            
                eye = eye[:,stim_on_time[0]-250:stim_off_time[0]-150]
                
                # get speed of eye movements
                dx = np.diff(eye[0,:])
                dy = np.diff(eye[1,:])
    
                #eye_speed = stats.zscore(np.hypot(dx,dy))
                eye_speed = stats.zscore(np.hypot(dx,dy))
    
                # find the saccades
                sacc_ix = sig.find_peaks(eye_speed,4,distance = 35)[0]
                
                # store some info about the saccades
                # how many saccades?
                sessiondf.loc[t,13] = len(sacc_ix)
                
                # loop over number of saccades
                for ix,s in enumerate(sacc_ix):
                    x_pos = np.mean(eye[0,s+1:s+10])
                    
                    if x_pos < -5:
                        sessiondf.loc[t,13+ix+1] = -1
                    elif x_pos > 5:
                        sessiondf.loc[t,13+ix+1] = 1
                    elif (x_pos < 5) & (x_pos > -5):
                            sessiondf.loc[t,13+ix+1] = 0
                            
                # get value associated with the first saccade
                if sessiondf.loc[t,14] == -1:
                   sessiondf.loc[t,19] = sessiondf.loc[t,7]
                else:
                   sessiondf.loc[t,19] = sessiondf.loc[t,8]
                   
                        
      
                # #plot eye movement for the trial
                # plt.figure()
                # colors = cm.magma(np.linspace(0,1,len(eye_speed)))
                
                # plt.subplot(1,2,1)
                # plt.scatter(eye[0,0:len(eye_speed)], eye[1,0:len(eye_speed)], color=colors)
                # plt.xlim([-15,15])
                # plt.ylim([-15,15])
                
                # plt.subplot(1,2,2)
                # xvals = np.arange(len(eye_speed))*2
                # plt.plot(xvals,eye[0,0:len(eye_speed)], label = 'eye x')
                # plt.plot(xvals,eye_speed, label = 'speed')
                # plt.scatter(xvals, np.ones(len(eye_speed))*-12, color=colors)
    
                # plt.ylim([-15,15])
                # plt.plot([sacc_ix*2,sacc_ix*2],plt.ylim(),color='black',linewidth=1)
                # plt.legend()
                # plt.show()
        
                xx=[]
            

        # add these new data to the table
        alldata =alldata.append(sessiondf)

        # check whether to save the data
        if save_dir is not None:

            sessiondf.rename({1:'fname', 2:'tnum',3:'use',4:'rule',5:'rtype',6:'forced',7:'lval',
                    8:'rval',9:'chosenval',10:'pickedbest',11:'rt',12:'side', 
                    13: 'n_saccs', 14: 'sacc1', 15: 'sacc2', 16: 'sacc3',
                    17: 'sacc4', 18: 'sacc5', 19: 'sacc1_val'}, axis= 1,inplace=True)

            save_name = save_dir + fnames[i][0:-3] + '_bhv.h5'

            # initialize an h5 file
            sessiondf.to_hdf(save_name, key='trialinfo', mode='w')

    # add some labels to the dataframe
    alldata.rename({1:'fname', 2:'tnum',3:'use',4:'rule',5:'rtype',6:'forced',7:'lval',
                    8:'rval',9:'chosenval',10:'pickedbest',11:'rt',12:'side', 
                    13: 'n_saccs', 14: 'sacc1', 15: 'sacc2', 16: 'sacc3',
                    17: 'sacc4', 18: 'sacc5', 19: 'sacc1_val'}, axis= 1,inplace=True)
    
    alldata = alldata.reset_index(drop=True)
    trials2drop = alldata.use == 0 # drop unusable trials
    alldata = alldata.drop(alldata[trials2drop].index)
    alldata = alldata.drop(labels='use',axis='columns')
    alldata = alldata.reset_index(drop=True)    
    
    xx=[]
    return alldata
# END of getMLH5



def assesschoice(alldata):
    
    # first, look at choice accuracy over sessions
        
    
    fnames = alldata.fname.unique()
    
    smeans = pd.DataFrame()
    
    
    for s in range(len(fnames)):
        fix = alldata.fname == fnames[s]
        
        smeans.at[s,'rule1mean'] = alldata.pickedbest[fix & (alldata.rule==1)].mean()
        smeans.at[s,'rule1sem'] = alldata.pickedbest[fix & (alldata.rule==1)].sem()
        smeans.at[s,'rule2mean'] = alldata.pickedbest[fix & (alldata.rule==2)].mean()
        smeans.at[s,'rule2sem'] = alldata.pickedbest[fix & (alldata.rule==2)].sem()

        
    # look at likelihood of choosing left option as fxn of offer diffs
    alldata['pickedleft'] = alldata['side'] == -1
    alldata['pickedleft']=alldata['pickedleft'].astype(int)
    alldata['offerdiff'] = alldata.lval - alldata.rval
    alldata["maxval"]    = np.maximum(alldata.lval,alldata.rval)    
        
    pchooseleft = pd.DataFrame()
    pchooseleft['rule1mean'] = alldata[alldata.rule==1].groupby('offerdiff').mean().pickedleft  
    pchooseleft['rule1sem'] = alldata[alldata.rule==1].groupby('offerdiff').sem().pickedleft  

    pchooseleft['rule2mean'] = alldata[alldata.rule==2].groupby('offerdiff').mean().pickedleft  
    pchooseleft['rule2sem'] = alldata[alldata.rule==2].groupby('offerdiff').sem().pickedleft    
    
    rt_x_cond = pd.DataFrame()
    rt_x_cond['rule1mean'] = alldata[alldata.rule==1].groupby('offerdiff').mean().rt  
    rt_x_cond['rule1sem'] = alldata[alldata.rule==1].groupby('offerdiff').sem().rt  

    rt_x_cond['rule2mean'] = alldata[alldata.rule==2].groupby('offerdiff').mean().rt  
    rt_x_cond['rule2sem'] = alldata[alldata.rule==2].groupby('offerdiff').sem().rt 
     
    valdiffs = alldata['offerdiff'].dropna().unique()
    valdiffs.sort()    
           
    s_xvals =  np.arange(0,len(fnames)); 
         
    fig, axs = plt.subplots()
    ax1 = plt.subplot2grid((2,2),(0,0))
    ax2 = plt.subplot2grid((2,2),(0,1))
    ax3 = plt.subplot2grid((2,2),(1,0))
    ax4 = plt.subplot2grid((2,2),(1,1))
    plt.tight_layout(pad = 2)
    
    
    # plot the session means
    ax1.errorbar(s_xvals,smeans['rule1mean'],smeans['rule1sem'],
                    color='tab:red',marker='o')
    ax1.errorbar(s_xvals,smeans['rule2mean'],smeans['rule2sem'],
                    color='tab:blue',marker='o')
    
    ax1.legend(['rule 1','rule 2'],loc='lower right')
    
    ax1.plot(s_xvals,np.ones(len(fnames))*.5,linestyle='--',color='gray')

    ax1.set_ylim(0,1)
    ax1.set_xticks(s_xvals)
    ax1.set_xticklabels(ax1.get_xticks()+1)
    ax1.set_xlabel('Session Num')
    ax1.set_ylabel('Percent Correct')
  
    

    # plot p(Choose Left)
    ax2.errorbar(valdiffs,pchooseleft['rule1mean'],pchooseleft['rule1sem'],
                    color = 'tab:red',marker='o')
    ax2.errorbar(valdiffs,pchooseleft['rule2mean'],pchooseleft['rule2sem'],
                    color = 'tab:blue',marker='o')
    
    ax2.plot(valdiffs,np.ones(len(valdiffs))*.5,linestyle='--',color='gray')
    ax2.set_xticks(np.arange(np.min(valdiffs),np.max(valdiffs)+1))
    
    ax2.set_xlabel('Lval - Rval')
    ax2.set_ylabel('p(Choose Left)')
    
    # let's do a logistic regression to summarize the offer diff data
    eqn = 'pickedleft ~ offerdiff*C(rule)'
    logitmodel = smf.glm(formula=eqn, data=alldata, family=sm.families.Binomial()).fit()
    # exctract the results
    modelresults = pd.DataFrame()
    modelresults['factor'] = pd.DataFrame(logitmodel.params.index)
    modelresults['coeff'] = logitmodel.params.values # this keeps only the value and not the index
    modelresults['pval'] = logitmodel.pvalues.values

    
    # break out the likelihood of choosing the better option by choice condition
    vals = alldata.rval.dropna().unique()
    vals.sort()
    
    choicegrid = np.empty((len(vals),len(vals)))
    choicegrid[:] = np.nan
    
    for val1 in range(len(vals)):
        for val2 in range(len(vals)):
            
            v1_ix = (alldata.rval == val1+1) | (alldata.lval == val1+1)
            v2_ix = (alldata.rval == val2+1) | (alldata.lval == val2+1)
            
            v1_ix = (alldata.rval == val1+1)
            v2_ix = (alldata.lval == val2+1)
            
            
            if val1 != val2:
                choicegrid[val1,val2] = np.nanmean(alldata.pickedbest[v1_ix&v2_ix])
            
    #choicegrid = np.triu(choicegrid)
    #choicegrid[choicegrid==0]=np.nan        
    gridim = ax3.imshow(choicegrid,vmin=0,vmax=1,origin = 'lower',cmap='seismic')
    ax3.set_xticks(np.arange(len(vals)))
    ax3.set_yticks(np.arange(len(vals)))
    fig.colorbar(gridim,ax =ax3)
    ax3.set_xlabel('val 1')
    ax3.set_ylabel('val 2')
    ax3.set_xticklabels(ax3.get_xticks()+1)
    ax3.set_yticklabels(ax3.get_xticks()+1)
    

    # plot choice RTs as function of value difference
    ax4.errorbar(valdiffs,rt_x_cond['rule1mean'],rt_x_cond['rule1sem'],
                        color = 'tab:red',marker='o')
    ax4.errorbar(valdiffs,rt_x_cond['rule2mean'],rt_x_cond['rule2sem'],
                        color = 'tab:blue',marker='o')
        
    ax4.set_xticks(np.arange(np.min(valdiffs),np.max(valdiffs)+1))
        
    ax4.set_xlabel('Lval - Rval')
    ax4.set_ylabel('rt (ms)')

    xx=[]
    return modelresults,choicegrid



def choice_or_rt_by_val_v01(alldata):
    '''
    A function to characterize self-control. 
    It looks at choices, RTs, and saccades as a function of an option
    value being present in a trial. 
    '''

    vals = alldata.rval.dropna().unique()
    vals.sort() 
    
    choice_means = np.empty((len(vals),2))
    choice_sems  = np.empty((len(vals),2)) 
    
    rt_means = np.empty((len(vals),2))
    rt_sems  = np.empty((len(vals),2)) 
    
    n_sacc_means = np.empty((len(vals),2))
    n_sacc_sems  = np.empty((len(vals),2)) 
    
    sacc1_means = np.empty((len(vals),2))
    sacc1_sems  = np.empty((len(vals),2)) 
    
    for v in range(len(vals)):

        val_present_ix = (alldata.lval == vals[v]) | (alldata.rval == vals[v])
        
        val_trial_data = alldata[val_present_ix]
        val_rule1_ix = val_trial_data.rule ==1

        chosen_ix = val_trial_data.chosenval == vals[v]
        picked4   = val_trial_data.chosenval == 4
        forced_ix = val_trial_data.forced ==1

        choice_means[v,0] = chosen_ix[val_rule1_ix].mean()
        choice_means[v,1] = chosen_ix[~val_rule1_ix].mean()
        choice_sems[v,0]  = chosen_ix[val_rule1_ix].sem()
        choice_sems[v,1]  = chosen_ix[~val_rule1_ix].sem()
        
        rt_means[v,0]     = val_trial_data.rt[val_rule1_ix].mean()
        rt_means[v,1]     = val_trial_data.rt[~val_rule1_ix].mean()
        rt_sems[v,0]      = val_trial_data.rt[val_rule1_ix].sem()
        rt_sems[v,1]      = val_trial_data.rt[~val_rule1_ix].sem()
        
        n_sacc_means[v,0] = val_trial_data.n_saccs[val_rule1_ix].mean()
        n_sacc_means[v,1] = val_trial_data.n_saccs[~val_rule1_ix].mean()
        n_sacc_sems[v,0]  = val_trial_data.n_saccs[val_rule1_ix].sem()
        n_sacc_sems[v,1]  = val_trial_data.n_saccs[~val_rule1_ix].sem()
        
        sacc1_means[v,0]  = val_trial_data.sacc1_val[val_rule1_ix].mean()
        sacc1_means[v,1]  = val_trial_data.sacc1_val[~val_rule1_ix].mean()
        sacc1_sems[v,0]   = val_trial_data.sacc1_val[val_rule1_ix].sem()
        sacc1_sems[v,1]   = val_trial_data.sacc1_val[~val_rule1_ix].sem()
        
        
        
    # create figure and axes 
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(6, 6))
    plt.tight_layout(pad = 2)
    
    ax1.errorbar(vals,choice_means[:,0],choice_sems[:,0],color = 'tab:red',marker='o')
    ax1.errorbar(vals,choice_means[:,1],choice_sems[:,1],color = 'tab:blue',marker='o')
    ax1.set_xticks(ticks=vals)
    ax1.set_xlabel('Option Value in Trial')
    ax1.set_ylabel('p(Choose Option)')    
    ax1.legend(['rule 1','rule 2'])
       
    ax2.errorbar(vals,rt_means[:,0],rt_sems[:,0],color = 'tab:red',marker='o')
    ax2.errorbar(vals,rt_means[:,1],rt_sems[:,1],color = 'tab:blue',marker='o')
    ax2.set_xticks(ticks=vals)
    ax2.set_xlabel('Option Value in Trial')
    ax2.set_ylabel('rt (ms)')    
    
    ax3.errorbar(vals,n_sacc_means[:,0],n_sacc_sems[:,0],color = 'tab:red',marker='o')
    ax3.errorbar(vals,n_sacc_means[:,1],n_sacc_sems[:,1],color = 'tab:blue',marker='o')
    ax3.set_xticks(ticks=vals)
    ax3.set_xlabel('Option Value in Trial')
    ax3.set_ylabel('# saccades') 
    
    ax4.errorbar(vals,sacc1_means[:,0],sacc1_sems[:,0],color = 'tab:red',marker='o')
    ax4.errorbar(vals,sacc1_means[:,1],sacc1_sems[:,1],color = 'tab:blue',marker='o')
    ax4.set_xticks(ticks=vals)
    ax4.set_xlabel('Option Value in Trial')
    ax4.set_ylabel('Value of First Saccade') 
    
    plt.show()
    
    xx=[]
    return
# END of choice_or_rt_by_val_v01

def choice_or_rt_by_val_v02(alldata):
    '''
    The difference v01 and v02 is that v02
    1. collapses across rules
    2. asks about the likelihood of the first saccade being towards option 1
    '''

    vals = alldata.rval.dropna().unique()
    vals.sort() 
    
    choice_means = np.empty((len(vals),1))
    choice_sems  = np.empty((len(vals),1)) 
    
    rt_means = np.empty((len(vals),1))
    rt_sems  = np.empty((len(vals),1))     
    
    sacc1_means = np.empty((len(vals),1))
    sacc1_sems  = np.empty((len(vals),1)) 
    
    for v in range(len(vals)):

        val_present_ix = (alldata.lval == vals[v]) | (alldata.rval == vals[v])
        
        val_trial_data = alldata[val_present_ix]

        chosen_ix = val_trial_data.chosenval == vals[v]
        picked4_ix   = val_trial_data.chosenval == 4
        forced_ix = val_trial_data.forced ==1
        firstsacc1 = val_trial_data.sacc1_val ==1

        choice_means[v,0] = chosen_ix.mean()
        choice_sems[v,0]  = chosen_ix.sem()
 
        rt_means[v,0]     = val_trial_data.rt[picked4_ix].mean()
        rt_sems[v,0]      = val_trial_data.rt[picked4_ix].sem()
               
        sacc1_means[v,0]  = firstsacc1[picked4_ix].mean()
        sacc1_sems[v,0]   = firstsacc1[picked4_ix].sem()
        
        
        
    # create figure and axes 
    fig, (ax1, ax2,ax3) = plt.subplots(1,3,figsize=(9, 3))
    plt.tight_layout(pad = 2)
    
    ax1.errorbar(vals,choice_means[:,0],choice_sems[:,0],color = 'black',marker='o')
    ax1.set_xticks(ticks=vals)
    ax1.set_xlabel('Option Value in Trial')
    ax1.set_ylabel('p(Choose Option)')    
       
    ax2.errorbar(vals,rt_means[:,0],rt_sems[:,0],color = 'black',marker='o')
    ax2.set_xticks(ticks=vals)
    ax2.set_xlabel('Option Value in Trial')
    ax2.set_ylabel('rt (ms)')    
    
    ax3.errorbar(vals,sacc1_means[:,0],sacc1_sems[:,0],color = 'black',marker='o')
    ax3.set_xticks(ticks=vals)
    ax3.set_xlabel('Option Value in Trial')
    ax3.set_ylabel('p(Sacc1 == 1)') 

    plt.show()
    
    xx=[]
    return
# END of choice_or_rt_by_val_v02


def compare_choice_rt_sacc1_x_conditions(normaldata,overrep_data):
    '''
    Compares choice, RT, and saccade patterns across conditions where the rules
    were 50/50 and 66/33
    '''
 

    vals = normaldata.rval.dropna().unique()
    vals.sort() 
    
    choice_means = np.empty((len(vals),3))
    choice_sems  = np.empty((len(vals),3)) 
    
    rt_means = np.empty((len(vals),3))
    rt_sems  = np.empty((len(vals),3))     
    
    sacc1_means = np.empty((len(vals),3))
    sacc1_sems  = np.empty((len(vals),3)) 
    
    all_rt_sacc_df = pd.DataFrame()
    all_choice_df = pd.DataFrame()
    
    for v in range(len(vals)):
        
        val_data = pd.DataFrame()

        norm_val_ix = (normaldata.lval == vals[v]) | (normaldata.rval == vals[v])
        overrep_val_ix = (overrep_data.lval == vals[v]) | (overrep_data.rval == vals[v])

        norm_val_trial_data = normaldata[norm_val_ix]
        overrep_val_trial_data = overrep_data[overrep_val_ix]

        norm_chosen_ix = norm_val_trial_data.chosenval == vals[v]
        norm_picked4_ix   = norm_val_trial_data.chosenval == 4
        norm_firstsacc1 = norm_val_trial_data.sacc1_val ==1
        
        overrep_chosen_ix = overrep_val_trial_data.chosenval == vals[v]
        overrep_picked4_ix   = overrep_val_trial_data.chosenval == 4
        overrep_firstsacc1 = overrep_val_trial_data.sacc1_val ==1
        overrep_rule1_ix = overrep_val_trial_data.rule == 1
        

        choice_means[v,0] = norm_chosen_ix.mean()
        choice_sems[v,0]  = norm_chosen_ix.sem()
        choice_means[v,1] = overrep_chosen_ix[overrep_rule1_ix].mean()
        choice_sems[v,1]  = overrep_chosen_ix[overrep_rule1_ix].sem()
        choice_means[v,2] = overrep_chosen_ix[~overrep_rule1_ix].mean()
        choice_sems[v,2]  = overrep_chosen_ix[~overrep_rule1_ix].sem()
 
        rt_means[v,0]     = norm_val_trial_data.rt[norm_picked4_ix].mean()
        rt_sems[v,0]      = norm_val_trial_data.rt[norm_picked4_ix].sem()
        rt_means[v,1]     = overrep_val_trial_data.rt[overrep_picked4_ix & overrep_rule1_ix].mean()
        rt_sems[v,1]      = overrep_val_trial_data.rt[overrep_picked4_ix & overrep_rule1_ix].sem()
        rt_means[v,2]     = overrep_val_trial_data.rt[overrep_picked4_ix & ~overrep_rule1_ix].mean()
        rt_sems[v,2]      = overrep_val_trial_data.rt[overrep_picked4_ix & ~overrep_rule1_ix].sem()
               
        sacc1_means[v,0]  = norm_firstsacc1[norm_picked4_ix].mean()
        sacc1_sems[v,0]   = norm_firstsacc1[norm_picked4_ix].sem()
        sacc1_means[v,1]  = overrep_firstsacc1[overrep_picked4_ix & overrep_rule1_ix].mean()
        sacc1_sems[v,1]   = overrep_firstsacc1[overrep_picked4_ix & overrep_rule1_ix].sem()
        sacc1_means[v,2]  = overrep_firstsacc1[overrep_picked4_ix & ~overrep_rule1_ix].mean()
        sacc1_sems[v,2]   = overrep_firstsacc1[overrep_picked4_ix & ~overrep_rule1_ix].sem()
        
        # build arrays for aggregation 
        # for choice data        
        choice_data = np.concatenate([norm_chosen_ix, 
                                      overrep_chosen_ix[overrep_rule1_ix],
                                      overrep_chosen_ix[~overrep_rule1_ix] ])
        
        choice_rules = np.concatenate([np.ones(len(norm_chosen_ix)),
                                 np.ones(sum(overrep_rule1_ix))*2,
                                 np.ones(sum(~overrep_rule1_ix))*3])
        
        choice_vals = np.concatenate([np.ones(len(norm_chosen_ix))*v,
                                      np.ones(sum(overrep_rule1_ix))*v,
                                      np.ones(sum(~overrep_rule1_ix))*v])
           
        # for rt data
        rt_data = np.concatenate([norm_val_trial_data.rt[norm_picked4_ix],
                                  overrep_val_trial_data.rt[overrep_picked4_ix & overrep_rule1_ix],
                                  overrep_val_trial_data.rt[overrep_picked4_ix & ~overrep_rule1_ix]])
        
        # for sacc1 data
        sacc1_data = np.concatenate([norm_firstsacc1[norm_picked4_ix],
                                     overrep_firstsacc1[overrep_picked4_ix & overrep_rule1_ix],
                                     overrep_firstsacc1[overrep_picked4_ix & ~overrep_rule1_ix]])
        
        # sacc and rt val and rule factors
        rt_sacc_rules = np.concatenate([np.ones(sum(norm_picked4_ix)),
                                        np.ones(sum(overrep_picked4_ix & overrep_rule1_ix))*2,
                                        np.ones(sum(overrep_picked4_ix & ~overrep_rule1_ix))*3])
        
        rt_sacc_vals = np.concatenate([np.ones(sum(norm_picked4_ix))*v,
                                        np.ones(sum(overrep_picked4_ix & overrep_rule1_ix))*v,
                                        np.ones(sum(overrep_picked4_ix & ~overrep_rule1_ix))*v])
        
        choice_df = pd.DataFrame()
        choice_df['choice'] = choice_data.astype(int)
        choice_df['rule']   = choice_rules 
        choice_df['val'] = (choice_vals+1).astype(int)
       
        # aggregate choice data into larger dataframe
        all_choice_df = all_choice_df.append(choice_df)
        
        sacc_rt_df = pd.DataFrame()
        sacc_rt_df['rt'] = rt_data 
        sacc_rt_df['sacc1'] = sacc1_data.astype(int)
        sacc_rt_df['rule'] = rt_sacc_rules 
        sacc_rt_df['val']  = (rt_sacc_vals+1).astype(int)
        
        # aggregate sacc_rt data into larger dataframe
        all_rt_sacc_df = all_rt_sacc_df.append(sacc_rt_df)
        
        xx=[]
        
        
        
    # create figure and axes 
    fig, (ax1, ax2,ax3) = plt.subplots(1,3,figsize=(9, 3), dpi=300)
    plt.tight_layout(pad = 2)
    
    ax1.errorbar(vals,choice_means[:,0],choice_sems[:,0],color = 'black',marker='o')
    ax1.errorbar(vals,choice_means[:,1],choice_sems[:,1],color = 'tab:red',marker='o')
    ax1.errorbar(vals,choice_means[:,2],choice_sems[:,2],color = 'tab:blue',marker='o')
    ax1.set_xticks(ticks=vals)
    ax1.set_xlabel('Option Value in Trial')
    ax1.set_ylabel('p(Choose Option)')    
       
    ax2.errorbar(vals,rt_means[:,0],rt_sems[:,0],color = 'black',marker='o')
    ax2.errorbar(vals,rt_means[:,1],rt_sems[:,1],color = 'tab:red',marker='o')
    ax2.errorbar(vals,rt_means[:,2],rt_sems[:,2],color = 'tab:blue',marker='o')
    ax2.set_xticks(ticks=vals)
    ax2.set_xlabel('Option Value in Trial')
    ax2.set_ylabel('rt (ms)')    
    
    ax3.errorbar(vals,sacc1_means[:,0],sacc1_sems[:,0],color = 'black',marker='o')
    ax3.errorbar(vals,sacc1_means[:,1],sacc1_sems[:,1],color = 'tab:red',marker='o')
    ax3.errorbar(vals,sacc1_means[:,2],sacc1_sems[:,2],color = 'tab:blue',marker='o')
    
    ax3.set_xticks(ticks=vals)
    ax3.set_xlabel('Option Value in Trial')
    ax3.set_ylabel('p(Sacc1 == 1)') 

    plt.show()
    #fig.savefig('ChoicexRTxSacc_v01.svg')
    
    
    # do stats
    choice_stats = pg.anova(data = all_choice_df,
                            dv = 'choice',
                            between = ['rule','val'])

    
    rt_stats     = pg.anova(data = all_rt_sacc_df,
                            dv = 'rt',
                            between = ['rule','val'])
    
    rt_posthoc = all_rt_sacc_df.pairwise_ttests(dv='rt', between=['val','rule'] )

    
    # since he couldn't look at 1 when it wasn't present, only look for those trials
    sacc1_keep_trials = (all_rt_sacc_df['val']==1) | (all_rt_sacc_df['val']==4)
    keep_sacc1_data = all_rt_sacc_df[sacc1_keep_trials]
    
    sacc1_stats  = pg.anova(data = keep_sacc1_data,
                            dv = 'sacc1',
                            between = ['rule','val'])
    
    # do some mulitple comparisions     
    sacc1_posthoc = all_rt_sacc_df.pairwise_ttests(dv='sacc1', between=['val','rule'] )
    
    
    xx=[]
    return
# END of compare_choice_rt_sacc1_x_conditions




def plot_rt_by_value(alldata):
        
    
    vals = alldata.rval.dropna().unique()
    vals.sort() 
    
    rt_means = np.empty((len(vals),2))
    rt_sems  = np.empty((len(vals),2)) 
    
    for v in range(len(vals)):

        val_present_ix = (alldata.lval == vals[v]) | (alldata.rval == vals[v])
        
        val_trial_data = alldata[val_present_ix]

        chosen_ix = val_trial_data.chosenval == vals[v]
        forced_ix = val_trial_data.forced ==1
        hit_ix    = val_trial_data.pickedbest ==1
        
        rt_means[v,0]     = val_trial_data.rt[hit_ix & chosen_ix & ~forced_ix].mean()
        rt_means[v,1]     = val_trial_data.rt[~hit_ix & chosen_ix & ~forced_ix].mean()
        rt_sems[v,0]      = val_trial_data.rt[hit_ix & chosen_ix & ~forced_ix].sem()
        rt_sems[v,1]      = val_trial_data.rt[~hit_ix & chosen_ix & ~forced_ix].sem()
        
    plt.figure    
    plt.errorbar(vals,rt_means[:,0],rt_sems[:,0],color = 'tab:blue',marker='o')
    plt.errorbar(vals,rt_means[:,1],rt_sems[:,1],color = 'tab:gray',marker='o')
    plt.xticks(ticks=vals)
    plt.xlabel('Chosen Option Value')
    plt.ylabel('rt (ms)')    
    plt.legend(['Correct','Error'])
    xx=[]
    
    return
# END of plot_rt_by_value

def check_switch_cost(alldata):
        
    fnames = alldata.fname.unique()
    

    switched_to_rule = np.array([])
    prior_acc = np.array([])
    switch_acc = np.array([])
    after_acc = np.array([])
    
    # loop over each session  
    for f in fnames:
        
        f_ix = alldata['fname']==f
        fdata = alldata.loc[f_ix,:]
        fdata = fdata.reset_index(drop=True)
        
        # find sequences of 1s and 2s
        lens,pos,ids = ut.find_sequences(fdata.rule)
        
        # find where at least 4 trials of the same rule occurred in a row
        valid_seqs = lens > 3
        
        switch_trials = pos[valid_seqs] + lens[valid_seqs]
        switch_trials = switch_trials[switch_trials+1<len(fdata)]
        
        switched_to_rule = np.append(switched_to_rule,fdata['rule'][switch_trials])
        
        prior_acc = np.append(prior_acc,fdata['pickedbest'][switch_trials-1])
        switch_acc = np.append(switch_acc,fdata['pickedbest'][switch_trials] )
        after_acc = np.append(after_acc,fdata['pickedbest'][switch_trials+1] )
        
    
    # let's do some plotting and stats
    x_data = pd.DataFrame()
    x_data['new_rule'] = np.concatenate([switched_to_rule,switched_to_rule])
    x_data['before_after'] = np.concatenate([np.ones(len(prior_acc))*-1, 
                                             np.ones(len(prior_acc))])
    x_data['acc'] = np.concatenate([prior_acc, switch_acc])
    x_data['switch_id'] = np.concatenate([np.arange(len(prior_acc)), 
                                          np.arange(len(prior_acc))])
     
    rule1_ix = x_data.new_rule ==1
    prior_ix = x_data.before_after == -1
    
    # summarize data for plotting
    rule_1_acc_means = np.array([x_data.acc[rule1_ix & prior_ix].mean(),
                                x_data.acc[rule1_ix & ~prior_ix].mean()])
                                      
    rule_1_acc_sems = np.array([x_data.acc[rule1_ix & prior_ix].sem(),
                               x_data.acc[rule1_ix & ~prior_ix].sem()])
                                     
    rule_2_acc_means = np.array([x_data.acc[~rule1_ix & prior_ix].mean(),
                                x_data.acc[~rule1_ix & ~prior_ix].mean()])
                                          
    rule_2_acc_sems = np.array([x_data.acc[~rule1_ix & prior_ix].sem(),
                               x_data.acc[~rule1_ix & ~prior_ix].sem()])                                     
                                      
                    
    fig, (ax1) = plt.subplots(1,1,figsize=(3,4))
    
    
    ax1.errorbar([1,2],rule_1_acc_means,rule_1_acc_sems,color='tab:red', marker='o')
    ax1.errorbar([1,2],rule_2_acc_means,rule_1_acc_sems,color='tab:blue', marker='o')
    ax1.set_xticks([1,2])
    ax1.set_xticklabels(['trial before','switch'])
    ax1.set_ylabel('Choice Accuracy')
    ax1.legend(['2-->1','1-->2'])

    
    xx=[] 
    return   

        
    
    

























