function [TrainIX, leftoverIX,n_per_cond] = BalanceTrainSet_v02(Trials2Balance,Params2Balance, n_trials2use)
%------------------------------------------------------------------------------
% INPUTS
% Trials2Balance   - logical array of the trials you want to balance
% Params2Balance   - cell array where each element is a vector of categorical
%                    parameters to balance (e.g. choice value and side)
%                    *** each {column} of Params2Balance must have the same
%                    number of elements as Trials2Balance

% OUTPUTS
% TrainIX          - trial indices of a fully balanced training set
% leftoverIX       - trial indices of trials not included in TrainIX

% NOTES
% you could loop over this function to produce different subsets for
% differently partioned folds in a classifier analysis.
%------------------------------------------------------------------------------
% Thomas Elston
% telston@nurhopsi.org
% 17 Jan 2021

TrainIX=zeros(size(Trials2Balance));
leftoverIX=[];

NumParams2Balance = numel(Params2Balance);
Params2Balance = cell2mat(Params2Balance);

[CondParams,~,CondIXs] = unique(Params2Balance(Trials2Balance,:),'rows','sorted');
NumConds = numel(CondParams(:,1));

% find number of occurences for each condition
tbl = tabulate(CondIXs);
if n_trials2use == 0
    
    MinNumOfTrials2Keep = min(tbl(:,2));
    n_per_cond = MinNumOfTrials2Keep;

else
    
    MinNumOfTrials2Keep = n_trials2use;
    n_per_cond = n_trials2use;
    
end

% randomly select MinNumOfTrials2Keep from each condition
for c = 1:NumConds
    
    ThisCondTrialsIX = ismember(Params2Balance,CondParams(c,:),'rows') & Trials2Balance;
    % shuffle the indices and take the first MinNumOfTrials2Keep
    CondIndices = find(ThisCondTrialsIX);
    ShuffledIndices=CondIndices(randperm(numel(CondIndices)));
    
    Trials2Keep    = ShuffledIndices(1:MinNumOfTrials2Keep);
    LeftoverTrials = ShuffledIndices(MinNumOfTrials2Keep+1:end);

    TrainIX(Trials2Keep)=1;
    leftoverIX(LeftoverTrials)=1;
end % of looping through conditions

TrainIX = find(TrainIX);
leftoverIX = find(leftoverIX);


end % of function