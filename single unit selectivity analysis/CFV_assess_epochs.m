function [outdata] = CFV_assess_epochs(cue_FRs, choice_FRs, c_FR2, bhv, u_name, OFC_ix, s)

outdata = table; 

reg_tbl = bhv; 
reg_tbl.trial_num = (1:size(cue_FRs,1))'; 
reg_tbl.cue_FRs = cue_FRs; 
reg_tbl.choice_FRs = choice_FRs; 
reg_tbl.c_FR2 = c_FR2; 

% regress for state information at cue
cue_state_mdl = fitglm(reg_tbl, 'cue_FRs ~ state + state_type + trial_num + state:state_type');

% regress for state and value information at choice
choice_mdl = fitglm(reg_tbl, 'choice_FRs ~ state + state_type + chosenval + trial_num + state:chosenval');

% choice_mdl = fitglm(reg_tbl, 'choice_FRs ~ state + chosenval + trial_num + state:chosenval');

% fit regressions for chosen value for each context, separately
state_A_mdl = fitglm(reg_tbl(reg_tbl.state == 1,:), 'c_FR2 ~ trial_num + chosenval');
state_B_mdl = fitglm(reg_tbl(reg_tbl.state == -1,:), 'c_FR2 ~ trial_num + chosenval');

outdata.u_name = u_name; 
outdata.subject = s; 
outdata.ofc = OFC_ix; 
outdata.cue_state_b = cue_state_mdl.Coefficients.Estimate(2);
outdata.cue_state_pval = cue_state_mdl.Coefficients.pValue(2);
outdata.cue_type_b = cue_state_mdl.Coefficients.Estimate(3);
outdata.cue_type_pval = cue_state_mdl.Coefficients.pValue(3);
outdata.cue_state_x_type_b = cue_state_mdl.Coefficients.Estimate(5);
outdata.cue_state_x_type_pval = cue_state_mdl.Coefficients.pValue(5);


% outdata.choice_state_b = choice_mdl.Coefficients.Estimate(2);
% outdata.choice_state_pval = choice_mdl.Coefficients.pValue(2);
% outdata.choice_type_b = choice_mdl.Coefficients.Estimate(3);
% outdata.choice_type_pval = choice_mdl.Coefficients.pValue(3);
% outdata.choice_state_x_type_b = choice_mdl.Coefficients.Estimate(5);
% outdata.choice_state_x_type_pval = choice_mdl.Coefficients.pValue(5);

outdata.choice_state_b = choice_mdl.Coefficients.Estimate(2);
outdata.choice_state_pval = choice_mdl.Coefficients.pValue(2);
outdata.choice_val_b = choice_mdl.Coefficients.Estimate(4);
outdata.choice_val_pval = choice_mdl.Coefficients.pValue(4);
outdata.choice_stateval_b = choice_mdl.Coefficients.Estimate(6);
outdata.choice_stateval_pval = choice_mdl.Coefficients.pValue(6);

outdata.state_A_b = state_A_mdl.Coefficients.Estimate(2);
outdata.state_A_pval = state_A_mdl.Coefficients.pValue(2);
outdata.state_B_b = state_B_mdl.Coefficients.Estimate(2);
outdata.state_B_pval = state_B_mdl.Coefficients.pValue(2);


end % of function


