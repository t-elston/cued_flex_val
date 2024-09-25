function [p, pe2] = compute_PW2(dep_var, factors, factor_labels)

dep_var = zscore(dep_var);

[p ,anova_tbl, stats] = anovan(dep_var, factors, 'display','off', 'model','interaction');

N_trials = numel(dep_var);

% calculate partial omega squared
[w2, pe2, VarNames] = calculatePartialW2_v01(anova_tbl, N_trials);

end % of function

