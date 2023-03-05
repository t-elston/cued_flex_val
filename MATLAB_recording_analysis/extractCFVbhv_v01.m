function [bhv] = extractCFVbhv_v01(datadir)

bhv=table;

FileData= dir([datadir '*.bhv2']);
[nFiles,~] = size(FileData);

for f = 1:nFiles
    
    % load the data
    thisFileData = mlread([datadir FileData(f).name]);
    
    [~,trialsinsession] = size(thisFileData);
    
    sessiondata = table;
    
    % extract each trial's data
    for t = 1:trialsinsession
        
        sessiondata.fname{t}           = FileData(f).name;
        sessiondata.UseTrial(t)        = thisFileData(t).UserVars.UseTrial;
        sessiondata.rule(t)            = thisFileData(t).UserVars.rule;
        sessiondata.pickedbest(t)      = thisFileData(t).UserVars.PickedBestOpt;
        sessiondata.side(t)            = thisFileData(t).UserVars.SideChosen;
        sessiondata.pickedleft(t)      = thisFileData(t).UserVars.pickedleft;
        sessiondata.Lval(t)            = thisFileData(t).UserVars.OptionVals(1);
        sessiondata.Rval(t)            = thisFileData(t).UserVars.OptionVals(2);
        sessiondata.optdiff(t)         = sessiondata.Lval(t) - sessiondata.Rval(t);
        sessiondata.chosenval(t)       = thisFileData(t).UserVars.ChosenVal;
        sessiondata.rt(t)              = thisFileData(t).UserVars.RT;
        sessiondata.forcedchoice(t)    = thisFileData(t).UserVars.forced;

    end % of cycling through trials
    
    % remove trials without a response
    sessiondata(~sessiondata.UseTrial,:)=[];
    
    % drop the UseTrial columns
    sessiondata.UseTrial = [];
    
    % accumulate the data
     bhv = [bhv; sessiondata];
    
end % of cyling through files

end % of function