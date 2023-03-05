function [sig_factors, sig_R2] = parallel_sliding_GLM(firing_rates, bhv, thresh)

[n_trials, n_times] = size(firing_rates);

bhv.state(bhv.state == 2) = -1;

reg_tbl = table;
reg_tbl.state = bhv.state;
reg_tbl.val = bhv.chosenval;
reg_tbl.state_val = bhv.state.*bhv.chosenval;
reg_tbl.side = bhv.side; 
reg_tbl.prior_state = circshift(reg_tbl.state,1);
reg_tbl.prior_val = circshift(reg_tbl.val,1);
reg_tbl.im_id = bhv.chosen_im;
reg_tbl.valsum = bhv.Lval + bhv.Rval;
reg_tbl.valdiff = abs(bhv.Lval - bhv.Rval); 
reg_tbl.hi_lo = (bhv.chosenval > 2)*2 - 1;
reg_tbl.trial_num(1:numel(bhv.state)) = 1:numel(bhv.state);

formula = 'fr ~ val*state';


for t = 1:n_times
    
    
%     %             reg_tbl.fr = firing_rates(:,t);
%     %             reg_tbl.shuff_fr = shuffle(firing_rates(:,t));
%     %
%     %             mdl = fitglm(reg_tbl, formula);
%     
%     % extract the significant bins
%     sig_factors(1,t) = mdl.Coefficients.pValue(2) < thresh;
%     sig_factors(2,t) = mdl.Coefficients.pValue(3) < thresh;
%     sig_factors(3,t) = mdl.Coefficients.pValue(4) < thresh;

[state_b,~,~,~,state_stats] = regress(firing_rates(:,t),[ones(size(bhv.state)), bhv.state]);
[val_b,~,~,~,val_stats] = regress(firing_rates(:,t),[ones(size(bhv.state)), bhv.chosenval, reg_tbl.hi_lo]);
[stateval_b,~,~,~,stateval_stats] = regress(firing_rates(:,t),[ones(size(bhv.state)),bhv.state.*bhv.chosenval]);

sig_factors(1,t) = state_stats(3) < thresh;
sig_factors(2,t) = val_stats(3) < thresh;
sig_factors(3,t) = stateval_stats(3) < thresh;

sig_R2(1,t) = state_stats(1);
sig_R2(2,t) = val_stats(1);
sig_R2(3,t) = stateval_stats(1);

    
end % of looping over times


end