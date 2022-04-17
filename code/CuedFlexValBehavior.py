# -*- coding: utf-8 -*-
"""
CuedFlexValBehavior
Loads and parses H5 files generated from NIMH monkeylogic
@author: Thomas Elston
"""

import os
import pandas as pd
import numpy as np
import pdb
import h5py
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import statsmodels.api as sm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

import utils as ut


# this command opens the debugger
                  # leave this commented for first runs, then uncomment for
                  # interactive debugging


def getMLH5(datadir,run_debug):
    
    
    if run_debug:
        pdb.set_trace() 
    
    # set the path so python is pointed at the data folder
    os.chdir(datadir)
  
    
    fnames = os.listdir(datadir)
    
    alldata = pd.DataFrame()
    print("Loading file#: ",end=" ")
    for i in range(len(fnames)):
        
        print(str(i) ,end=" ")
                
        f = h5py.File(fnames[i],'r')
        
        ftrials  = list(f["ML"].keys())[1:-1]        
        sessiondf = pd.DataFrame()
      
        for t in range(len(ftrials)):
                    
            sessiondf.loc[t,1] = fnames[i]
            sessiondf.loc[t,2] = t
            sessiondf.loc[t,3] = f["ML"][ftrials[t]]["UserVars"]["UseTrial"][0]
            sessiondf.loc[t,4] = f["ML"][ftrials[t]]["UserVars"]["rule"][0]
            sessiondf.loc[t,5] = f["ML"][ftrials[t]]["UserVars"]["rtype"][0]
            sessiondf.loc[t,6] = f["ML"][ftrials[t]]["UserVars"]["forced"][0]
            
            if sessiondf.loc[t,6] == 0:
            
                sessiondf.loc[t,7] = f["ML"][ftrials[t]]["UserVars"]["OptionVals"][0]
                sessiondf.loc[t,8] = f["ML"][ftrials[t]]["UserVars"]["OptionVals"][1]
            else:
                if f["ML"][ftrials[t]]["UserVars"]["SideChosen"][0] == -1:
                    sessiondf.loc[t,7] = f["ML"][ftrials[t]]["UserVars"]["OptionVals"][0]
                    sessiondf.loc[t,8] = np.NaN
                else:
                    sessiondf.loc[t,8] = f["ML"][ftrials[t]]["UserVars"]["OptionVals"][0]
                    sessiondf.loc[t,7] = np.NaN
                    
    
            sessiondf.loc[t,9] = f["ML"][ftrials[t]]["UserVars"]["ChosenVal"][0]
            sessiondf.loc[t,10] = f["ML"][ftrials[t]]["UserVars"]["PickedBestOpt"][0]
            sessiondf.loc[t,11] = f["ML"][ftrials[t]]["UserVars"]["RT"][0]
            sessiondf.loc[t,12] = f["ML"][ftrials[t]]["UserVars"]["SideChosen"][0]            

            

        # add these new data to the table
        alldata =alldata.append(sessiondf)
        
    # add some labels to the dataframe
    alldata.rename({1:'fname', 2:'tnum',3:'use',4:'rule',5:'rtype',6:'forced',7:'lval',
                    8:'rval',9:'chosenval',10:'pickedbest',11:'rt',12:'side'}, axis= 1,inplace=True)
    
    trials2drop = alldata.use !=1  
    alldata = alldata.drop(alldata[trials2drop].index)
    alldata = alldata.drop(labels="use",axis='columns')
        
  
    return alldata



def assesschoice(alldata,run_debug):
    
    # first, look at choice accuracy over sessions
    if run_debug:
        pdb.set_trace() 
        
    
    fnames = alldata.fname.unique()
    
    smeans = pd.DataFrame()
    
    
    for s in range(len(fnames)):
        fix = alldata.fname == fnames[s]
        
        smeans.at[s,'rule1mean'] = alldata.pickedbest[fix & (alldata.rule==1)].mean()
        smeans.at[s,'rule1sem'] = alldata.pickedbest[fix & (alldata.rule==1)].sem()
        smeans.at[s,'rule2mean'] = alldata.pickedbest[fix & (alldata.rule==2)].mean()
        smeans.at[s,'rule2sem'] = alldata.pickedbest[fix & (alldata.rule==2)].sem()

        
    # look at likelihood of choosing left option as fxn of offer diffs
    alldata["pickedleft"] = alldata["side"] == -1
    alldata['pickedleft']=alldata['pickedleft'].astype(int)
    alldata["offerdiff"] = alldata.lval - alldata.rval
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
         
    #%% main plot
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



def choice_or_rt_by_val(alldata,debug_):
 #%% forced choice rt plot
 
    if debug_:
        pdb.set_trace() 
 

    vals = alldata.rval.dropna().unique()
    vals.sort() 
    
    choice_means = np.empty((len(vals),2))
    choice_sems  = np.empty((len(vals),2)) 
    
    rt_means = np.empty((len(vals),2))
    rt_sems  = np.empty((len(vals),2)) 
    
    for v in range(len(vals)):

        val_present_ix = (alldata.lval == vals[v]) | (alldata.rval == vals[v])
        
        val_trial_data = alldata[val_present_ix]
        val_rule1_ix = val_trial_data.rule ==1

        chosen_ix = val_trial_data.chosenval == vals[v]
        forced_ix = val_trial_data.forced ==1

        choice_means[v,0] = chosen_ix[val_rule1_ix].mean()
        choice_means[v,1] = chosen_ix[~val_rule1_ix].mean()
        choice_sems[v,0]  = chosen_ix[val_rule1_ix].sem()
        choice_sems[v,1]  = chosen_ix[~val_rule1_ix].sem()
        
        rt_means[v,0]     = val_trial_data.rt[val_rule1_ix].mean()
        rt_means[v,1]     = val_trial_data.rt[~val_rule1_ix].mean()
        rt_sems[v,0]      = val_trial_data.rt[val_rule1_ix].sem()
        rt_sems[v,1]      = val_trial_data.rt[~val_rule1_ix].sem()
        
        
        
    # create figure and axes 
    fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(6, 3))
    plt.tight_layout(pad = 2)
    
    ax1.errorbar(vals,choice_means[:,0],choice_sems[:,0],color = 'tab:red',marker='o')
    ax1.errorbar(vals,choice_means[:,1],choice_sems[:,1],color = 'tab:blue',marker='o')
    ax1.set_xticks(ticks=vals)
    ax1.set_xlabel('Option Value')
    ax1.set_ylabel('p(Choose Option)')    
    ax1.legend(['rule 1','rule 2'])
    
    
    ax2.errorbar(vals,rt_means[:,0],rt_sems[:,0],color = 'tab:red',marker='o')
    ax2.errorbar(vals,rt_means[:,1],rt_sems[:,1],color = 'tab:blue',marker='o')
    ax2.set_xticks(ticks=vals)
    ax2.set_xlabel('Option Value')
    ax2.set_ylabel('rt (ms)')    
    
    xx=[]
    return


























