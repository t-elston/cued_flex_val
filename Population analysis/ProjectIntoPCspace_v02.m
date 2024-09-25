function [projectedPCs] = ProjectIntoPCspace_v02(inspikes,mu,coeff,PC95)
projectedPCs=[];

   [nTrials,nTimeSteps,~] = size(inspikes);
    spikes4PCA = inspikes;
    term1 = spikes4PCA - repmat(mu,nTrials,1);
    term2 = repmat(coeff(:,1:PC95),1,1);   
    projectedPCs = pagemtimes(term1,term2); % NumTrials x pcs x NumTimeSteps 
    
    % rearrange
     projectedPCs=squeeze(permute(projectedPCs,[1,3,2])); 
    
end % of function