function [HPC_accs, OFC_accs] = pseudo_LDA_v01(train_labels, train_firing_rates, test_labels, test_firing_rates, OFC_ix, n_boots)

n_times = numel(train_firing_rates{1}(1,:,1));

HPC_accs = [];
OFC_accs = [];

for b = 1: n_boots
    
    HPC_train = train_firing_rates{b}(:,:,:);
    OFC_train = train_firing_rates{b}(:,:,OFC_ix ==1);
    
    rand_trials = shuffle(1:numel(HPC_train(:,1,1)));
    train_ix = rand_trials(1:round(.9*numel(HPC_train(:,1,1))));
    test_ix  = rand_trials(round(.9*numel(HPC_train(:,1,1))):end);
    
    for t = 1: n_times      
        
        HPC_pcs = getBestPCs_v01(squeeze(HPC_train(:,t,:)));
        OFC_pcs = getBestPCs_v01(squeeze(OFC_train(:,t,:)));

        HPC_LDA = fitcdiscr(HPC_pcs(train_ix,:),train_labels(train_ix,1));
        OFC_LDA = fitcdiscr(OFC_pcs(train_ix,:),train_labels(train_ix,1));
        
        HPC_pred = predict(HPC_LDA, HPC_pcs(test_ix,:));
        OFC_pred = predict(OFC_LDA, OFC_pcs(test_ix,:));
        
        HPC_accs(b, t) = mean(HPC_pred == train_labels(test_ix,1));
        OFC_accs(b, t) = mean(OFC_pred == train_labels(test_ix,1));

        
        
    end % of looping over timesteps
end % of looping over bootstraps

end % of function