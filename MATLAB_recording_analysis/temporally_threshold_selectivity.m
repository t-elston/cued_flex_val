function [new_sig_factors] = temporally_threshold_selectivity(sig_factors, thresh)

% go through each unit and find stretchs of selectivity that are at least
% thresh long

new_sig_factors=sig_factors;

[n_factors, n_times, n_units] = size(sig_factors);

for u = 1:n_units
    u_sig_factors = zeros(n_factors, n_times);
    for f = 1:n_factors
        [SelectivityStart, SelectivityLen, numEpochs] = ZeroOnesCount(sig_factors(f,:,u));
        
        if any(SelectivityLen >= thresh)
            % get index of when it started
           % window_start = SelectivityStart(find(SelectivityLen >= thresh));
            %win_len = SelectivityLen(find(SelectivityLen >= thresh));
            
            %u_sig_factors(f,window_start:window_start+win_len-1) = 1;
            
        else
            
            new_sig_factors(f,:,u) = 0;
            
            
        end
        
    end
    
    % aggregate the units 
    %new_sig_factors = cat(3,new_sig_factors, u_sig_factors);
end

end