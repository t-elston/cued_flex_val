function [cue_pvals, cue_betas, cue_sig, cue_aucs, choice_pvals, choice_betas] = screen_cue_and_choice_v01(firing_rates, bhv, ts)


% zscore the firing rates
fr_mean = nanmean(firing_rates, 'all');
fr_std  = nanstd(firing_rates,[], 'all');

firing_rates = (firing_rates - fr_mean) / fr_std;

cue_times_ix = ts >= -650 & ts <= -300;
choice_time_ix = ts >= 10 & ts <= 400;

cue_frs = nanmean(firing_rates(:,cue_times_ix),2);
choice_frs = nanmean(firing_rates(:,choice_time_ix),2);

cue_frs = cue_frs + abs(min(cue_frs));
choice_frs = choice_frs + abs(min(choice_frs));

bhv.Rval(isnan(bhv.Rval)) = 0;
bhv.Lval(isnan(bhv.Lval)) = 0;


bhv.state(bhv.state==2) = -1;

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
reg_tbl.state_image = bhv.state_type.*bhv.state;
reg_tbl.state_type = (bhv.state_type==1)*2 -1 ;


reg_tbl.cue_fr = cue_frs;
reg_tbl.choice_fr = choice_frs;

cue_formula = 'cue_fr ~ state*state_type + trial_num';
choice_formula = 'choice_fr ~ val*state + valdiff + trial_num + side + im_id';


cue_mdl = fitglm(reg_tbl, cue_formula);
choice_mdl = fitglm(reg_tbl, choice_formula);

[cue_p, cue_anova] = anovan(cue_frs,{reg_tbl.state, reg_tbl.state_type},'display','off','varnames',{'state','type'},'model','interaction'); 

% cue_pvals(1,1) = cue_mdl.Coefficients.pValue(2);
% cue_pvals(1,2) = cue_mdl.Coefficients.pValue(4);
% cue_pvals(1,3) = cue_mdl.Coefficients.pValue(5);
cue_pvals(1,:) = cue_p;

cue_betas(1,1) = cue_mdl.Coefficients.Estimate(2);
cue_betas(1,2) = cue_mdl.Coefficients.Estimate(4);
cue_betas(1,3) = cue_mdl.Coefficients.Estimate(5);



cue_sig = (cue_pvals(1,1) < .05) & (cue_pvals(1,2) > .05);


choice_pvals(1:2,1) = choice_mdl.Coefficients.pValue(2:3);
choice_betas(1:2,1) = choice_mdl.Coefficients.Estimate(2:3);
choice_pvals(3,1) = choice_mdl.Coefficients.pValue(8);
choice_betas(3,1) = choice_mdl.Coefficients.Estimate(8);



% AUROC when states match
[~,~,~,cue_aucs(1,1)] = perfcurve(reg_tbl.state,cue_frs,1);

% AUROC when types match
[~,~,~,cue_aucs(1,2)] = perfcurve(reg_tbl.state_type,cue_frs,1);

% AUROC when neither type nor state matches
s1_t1_frs = cue_frs(bhv.state_type==1 & bhv.state ==1);
s1_t2_frs = cue_frs(bhv.state_type==2 & bhv.state ==1);
s2_t1_frs = cue_frs(bhv.state_type==1 & bhv.state ==-1);
s2_t2_frs = cue_frs(bhv.state_type==2 & bhv.state ==-1); 

un_matched_frs = [s1_t1_frs ; s2_t2_frs ; s1_t2_frs ; s2_t1_frs];
un_match_lbls = [zeros(size(s1_t1_frs)) ; zeros(size(s2_t2_frs)) ; ones(size(s1_t2_frs)) ; ones(size(s2_t1_frs)) ];

[~,~,~,cue_aucs(1,3)] = perfcurve(un_match_lbls,un_matched_frs,1);




end % of function 



