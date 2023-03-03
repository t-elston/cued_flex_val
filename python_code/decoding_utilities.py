'''
decoding_utilities.py
@author: Thomas Elston
'''

import numpy as np
from scipy import stats

def get_labelled_posteriors(indata, labels):

    '''
    INPUTS:
    indata = posterior probabilites from a classifier with the shape
            n_trials x n_timesteps x n_classes
        
    labels = 1d array with len(n_trials) - these labels ought
            to correspond to class numbers (layers in indata)

    OUTPUT:
        labelled_posteriors = posterior probabilities associated with the
        classes in the labels input for each timestep and trial
    '''

    n_trials, n_times, n_classes = indata.shape
    class_lbls = np.unique(labels)
    class_lbls = class_lbls[~np.isnan(class_lbls)]

    # initialize output
    labelled_posteriors = np.zeros(shape = (n_trials, n_times))

    for ix, lbl in enumerate(class_lbls):
        
        # find trials where this label was chosen
        labelled_posteriors[labels == lbl,:] = indata[labels == lbl,:,int(ix)]
        
    return labelled_posteriors


def pull_balanced_train_set(trials2balance, params2balance):
    '''
    INPUTS:
    trials2balance   - ***logical array*** of the trials you want to balance
    params2balance   - ***list*** where each element is a vector of categorical
                        parameters to balance (e.g. choice value and side)
                        each element of params2balance must have the same
                        number of elements as trials2balance
    OUTPUTS:
    train_ix         - trial indices of a fully balanced training set
    leftover_ix      - trial indices of trials not included in train_ix

    NOTES:
    you could loop over this function to produce different subsets for
    differently partioned folds in a classifier analysis.
    '''

    # initialize output
    train_ix = np.zeros(shape = (len(trials2balance, )))
    leftover_ix =  np.zeros(shape = (len(trials2balance, )))

    # find out how many things to balance
    n_params2balance = len(params2balance)

    # convert params2balance to an array; .Transpose so rows are trials
    params2balance = np.array(params2balance).T.astype(int)

    # find the combos parameters and how many there are of each combo
    p_combos, p_counts = np.unique(params2balance[trials2balance,:], axis = 0, return_counts=True)
    
    # how many trials do we need to keep?
    n_to_keep = np.min(p_counts)

    # randomly select n_to_keep from each parameter
    for p in range(len(p_counts)):

        # find the indices where this current combo occurs
        this_param_ix = np.where((params2balance == p_combos[p,:]).all(axis = 1) & trials2balance)

        # convert from tuple to array
        this_param_ix = np.asarray(this_param_ix).flatten()

        # shuffle the array
        shuffled_indices = np.random.permutation(this_param_ix)

        # keep a number of the random trials
        trials2keep = shuffled_indices[0:n_to_keep]
        leftover_trials = shuffled_indices[n_to_keep: ]

        train_ix[trials2keep] = 1
        leftover_ix[leftover_trials] = 1
        
    # return train_ix as a boolean
    train_ix = train_ix == 1
    return train_ix


def random_prop_of_array(inarray, proportion):
    '''
    INPUTS
    inarray = logical/boolean array of indices to potentially use later
    proportion = how much of inarray should randomly be selected

    OUTPUT
    out_array = logical/boolean that's set as 'true' for a proportion of the 
                initial 'true' values in inarray
    '''

    out_array = np.zeros(shape = (len(inarray), ))

    # find where inarray is true and shuffle those indices
    shuffled_ixs = np.random.permutation(np.asarray(np.where(inarray)).flatten())

    # keep only a proportion of that array
    kept_ix = shuffled_ixs[0: round(len(shuffled_ixs)*proportion)]

    # fill in the kept indices
    out_array[kept_ix] = 1

    # make this a logical/boolean
    out_array = out_array > 0

    return out_array

def zscore_tensor(indata):
    '''
    INPUT
    indata = a tensor of data shaped n_trials x n_times x n_sensors/neurons

    OUTPUT
    outdata = a tensor where zscores have been computer for each sensor/neuron
    '''

    n_trials, n_times, n_units = indata.shape

    outdata = np.zeros_like(indata)

    outdata[:,:,np.arange(n_units)] = stats.zscore(indata[:,:,np.arange(n_units)], axis = None)

    return outdata

def ecdf(data):
    ''' Compute ECDF '''
    x = np.sort(data)
    n = x.size
    y = np.arange(1, n+1) / n
    return x,y











