function [z_data] = z_score_tensor_v01(indata)

% zscores along 3rd dim of indata

z_data = [];

for i = 1:numel(indata(1,1,:))
    
    i_data = indata(:,:,i);
    i_mean = nanmean(i_data,'all');
    i_std  = nanstd(i_data,[],'all');
    
    z_data(:,:,i) = (i_data - i_mean) / i_std;
    
end % of looping over layers of indata

end % of function