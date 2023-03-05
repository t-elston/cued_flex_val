function [bhv, ml_align_times] = extractCFVbhv_vRecorded(datadir, filename)

bhv=table;
ml_align_times=[];

warning('off','MATLAB:table:RowsAddedExistingVars');

% load the data
thisFileData = mlread([datadir filename]);

[~,trialsinsession] = size(thisFileData);

chosen_im = {};

% extract each trial's data
for t = 1:trialsinsession
    
    bhv.fname{t}           = filename(1:12);
    bhv.UseTrial(t)        = thisFileData(t).UserVars.UseTrial;
    bhv.state(t)           = thisFileData(t).UserVars.rule;
    bhv.state_type(t)      = thisFileData(t).UserVars.rtype;
    bhv.pickedbest(t)      = thisFileData(t).UserVars.PickedBestOpt;
    bhv.side(t)            = thisFileData(t).UserVars.SideChosen;
    bhv.pickedleft(t)      = thisFileData(t).UserVars.pickedleft;
    bhv.forcedchoice(t)    = thisFileData(t).UserVars.forced;
    bhv.condition(t)       = thisFileData(t).Condition;
    
    % was it a forced choice trial?
    if bhv.forcedchoice(t)
        bhv.Lval(t)            = NaN;
        bhv.Rval(t)            = NaN;
        bhv.optdiff(t)         = NaN;
        
        % was it forced left or right?
        if bhv.pickedleft(t) == 1
            bhv.Lval(t)= thisFileData(t).UserVars.OptionVals;
            chosen_im{t} = thisFileData(t).TaskObject.Attribute{1,3}{1,2};
        else
            bhv.Rval(t)= thisFileData(t).UserVars.OptionVals;
            chosen_im{t} = thisFileData(t).TaskObject.Attribute{1,2}{1,2};
        end % of parsing value for left/right forced trials
    else
        bhv.Lval(t) = thisFileData(t).UserVars.OptionVals(1);
        bhv.Rval(t) = thisFileData(t).UserVars.OptionVals(2);
        bhv.optdiff(t) = bhv.Lval(t) - bhv.Rval(t);
        
        
        if bhv.pickedleft(t) == 1
            chosen_im{t} = thisFileData(t).TaskObject.Attribute{1,3}{1,2};
        else
            chosen_im{t} = thisFileData(t).TaskObject.Attribute{1,2}{1,2};
        end
        
    end % of parsing choice on forced choice trials
    
    
    bhv.chosenval(t) = thisFileData(t).UserVars.ChosenVal;
    bhv.rt(t) = thisFileData(t).UserVars.RT; 
    
    % find out the absolute trial time of the ml align_event
    align_ix = max(find(thisFileData(t).BehavioralCodes.CodeNumbers == 3));
    try
    ml_align_times(t) = round((thisFileData(t).BehavioralCodes.CodeTimes(align_ix) + thisFileData(t).AbsoluteTrialStartTime));
    catch
        ml_align_times(t) = NaN;
    end
    
    
    
end % of cycling through trials

bhv.chosen_im = grp2idx(chosen_im);


end % of function