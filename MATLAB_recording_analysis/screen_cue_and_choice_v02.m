function [cue_pvals, cue_betas, cue_sig, choice_pvals, choice_betas] = screen_cue_and_choice_v02(firing_rates, bhv, ts)



cue_times_ix = ts >= -650 & ts <= -300;
choice_time_ix = ts >= 100 & ts <= 400;

cue_frs = nanmean(firing_rates(:,cue_times_ix),2);
choice_frs = nanmean(firing_rates(:,choice_time_ix),2);

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
reg_tbl.state_type = bhv.state_type;

reg_tbl.cue_fr = cue_frs;
reg_tbl.choice_fr = choice_frs;

cue_formula = 'cue_fr ~ state*state_type + trial_num';
choice_formula = 'choice_fr ~ val*state + valdiff + trial_num + side + im_id';

[cue_pvals, cue_anova] = anovan(cue_frs,{reg_tbl.state, reg_tbl.state_type},'model','interaction','display','off','varnames',{'state','type'});


cue_mdl = fitglm(reg_tbl, cue_formula);
choice_mdl = fitglm(reg_tbl, choice_formula);

cue_pvals(1,1) = cue_mdl.Coefficients.pValue(1);
cue_pvals(1,2) = cue_mdl.Coefficients.pValue(3);

cue_betas(1,1) = cue_mdl.Coefficients.Estimate(2);
[~,~,~,cue_betas(1,2)] = perfcurve(reg_tbl.state,cue_frs,1);
cue_betas(1,2) = cue_betas(1,2)*2 - 1;

cue_sig = (cue_pvals(1,1) < .05) & (cue_pvals(1,3) > .05);


choice_pvals(1:2,1) = choice_mdl.Coefficients.pValue(2:3);
choice_betas(1:2,1) = choice_mdl.Coefficients.Estimate(2:3);
choice_pvals(3,1) = choice_mdl.Coefficients.pValue(8);
choice_betas(3,1) = choice_mdl.Coefficients.Estimate(8);

end % of function 



