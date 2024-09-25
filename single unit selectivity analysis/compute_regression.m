function [p, beta] = compute_regression(dep_var, factor_table)

p=[];
beta=[];

factor_table.fr = zscore(dep_var);

mdl = fitlm(factor_table, 'fr ~ state*value');

p = mdl.Coefficients.pValue(2:end);
beta = mdl.Coefficients.Estimate(2:end);

end % of function