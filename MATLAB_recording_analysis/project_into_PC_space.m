function [projectedPCs] = project_into_PC_space(inspikes, mu,coeff)
projectedPCs=[];

   [nTrials,nTimeSteps,~] = size(inspikes);
    spikes4PCA = permute(inspikes,[1,3,2]);
    term1 = spikes4PCA - repmat(mu,nTrials,1,nTimeSteps);
    term2 = repmat(coeff,1,1,nTimeSteps);   
    projectedPCs = pagemtimes(term1,term2); % NumTrials x pcs x NumTimeSteps 
    
    % rearrange
      projectedPCs=(permute(projectedPCs,[1,3,2])); 
    
end % of function