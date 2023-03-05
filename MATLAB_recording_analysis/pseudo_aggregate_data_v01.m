function [train_labels, b_train_firing_rates, test_labels, b_test_firing_rates, OFC_ix] = pseudo_aggregate_data_v01(n_boots, datadir)

% find the files and their names
rec_files = dir([datadir '*.mat']);
n_files = numel(rec_files);

% load all firing rates and behavioral data precisely once 
b_train_firing_rates = {};
test_firing_rates = {};
train_labels={};
test_labels={};

for b = 1:n_boots
    
    boot_train_labels = [];
    boot_train_firing_rates = [];
    boot_test_firing_rates=[];
    
    OFC_ix = [];
    
    for f = 1:n_files
        
        % behavior
        bhv = load([datadir rec_files(f).name], 'bhv');
        bhv = bhv.bhv;
        
        bhv = bhv(bhv.UseTrial == 1,:);
        
        % firing_rates
        firing_rates = load([datadir rec_files(f).name], 'down_firing_rates');
        firing_rates = z_score_tensor_v01(firing_rates.down_firing_rates);
        firing_rates = firing_rates(bhv.UseTrial == 1,:,:);
        
        % sorting notes
        notes = load([datadir rec_files(f).name], 'sorting_notes');
        notes = notes.sorting_notes;
        
        % aggregate the
        % get train and test sets for each file
        [TrainIX, leftoverIX, n_per_cond] = BalanceTrainSet_v02(bhv.forcedchoice,{bhv.state, bhv.chosenval},10);
        [TestIX, leftoverIX, n_per_cond] = BalanceTrainSet_v02(bhv.forcedchoice ==0,{bhv.state, bhv.chosenval},10);

        
        % aggregate the labels and firing rates
        [train_labels{b}, train_firing_rates] = pseudo_aggregate_labels(TrainIX, bhv, firing_rates);
        [test_labels{b}, test_firing_rates] = pseudo_aggregate_labels(TestIX, bhv, firing_rates);

        
        boot_train_firing_rates = cat(3,boot_train_firing_rates, train_firing_rates);
        boot_test_firing_rates = cat(3,boot_test_firing_rates, test_firing_rates);
        
        OFC_ix = [OFC_ix; double(contains(notes.Target,'OFC'))];
        
    end % of looping over files
    
    b_train_firing_rates{b} = boot_train_firing_rates;
    b_test_firing_rates{b} = boot_test_firing_rates;
    
end % of looping over bootstraps




end % of function