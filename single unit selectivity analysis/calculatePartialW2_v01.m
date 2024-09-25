function [w2, pe2, VarNames] = calculatePartialW2_v01(anova_tbl,N_trials)
w2=[];

% def of partial omega-squared as per
% Maxwell, S. E., & Delaney, H. D. (2004). Designing experiments and analyzing data: A model comparison perspective (2nd ed.).
% Mahwah, New Jersey: Lawrence Erlbaum Associates

% and 

% Keren, G., & Lewis, C. (1979). Partial omega squared for ANOVA designs. Educational and Psychological Measurement, 39, 119–
% 128.

% w2 =   [DF_effect*(MS_effect - MS_error) ] 
%                      / 
%        [SS_effect + (N_total - DF_effect)*MS_error]
%     





% find the factor names between the 'source' heading and the 'error' footer

ANOVA_headers = anova_tbl(1,:);
ANOVA_SourceRows = anova_tbl(:,1);
ERROR_row        = find(strcmp(ANOVA_SourceRows,'Error'));

% find column with factor names
sourceHeader_location = find(strcmp(ANOVA_headers,'Source'));
SS_col = find(strcmp(ANOVA_headers,'Sum Sq.'));
DF_col = find(strcmp(ANOVA_headers,'d.f.'));
MS_col = find(strcmp(ANOVA_headers,'Mean Sq.'));


StartOfVarNames =  find(circshift(strcmp(anova_tbl(:,sourceHeader_location),'Source'),1));
EndOfVarNames =  find(circshift(strcmp(anova_tbl(:,sourceHeader_location),'Error'),-1));

VarNames = anova_tbl(StartOfVarNames:EndOfVarNames,sourceHeader_location);

% cycle through VarNames and calculate the partial-omega-squared

for v_ix = 1:numel(VarNames)
    thisVar = VarNames(v_ix);
    
    % find the row for this var
    VarRow = find(strcmp(ANOVA_SourceRows,thisVar{1}));
    
    % assemble the pieces for the partial omega squared
    Var_DF   = anova_tbl{VarRow,DF_col};
    Var_MS   = anova_tbl{VarRow,MS_col};
    Var_SS   = anova_tbl{VarRow,SS_col};
    Var_MS   = anova_tbl{VarRow,MS_col};
    Error_SS = anova_tbl{ERROR_row,SS_col};
    
    N_total  = N_trials;    
    MS_error = anova_tbl{ERROR_row,MS_col};
       
    partialOmegaSquared = (Var_DF*(Var_MS-MS_error)) / (Var_SS + ((N_total-Var_DF)*MS_error));
    partialEtaSquared   = Var_SS / (Var_SS + Error_SS);
    
    if partialOmegaSquared < 0
        partialOmegaSquared = 0;
        partialEtaSquared   = 0;
    end


    
    w2(v_ix,1) = partialOmegaSquared;
    pe2(v_ix,1) = partialEtaSquared;
    
end % of calculating w2 for each var


return