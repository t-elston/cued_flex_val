function [outIX] = get_random_proportion_v02(indata,proportion)

% this function works on one 1-D vector at a time
% Thomas Elston 
% telston@nurhopsi.org
% May 5 2021

outIX=[];

valid_ix = find(indata);

while isempty(outIX)
    
    random_indices = randperm(numel(valid_ix));
    keepIXs = random_indices(1:ceil(proportion*numel(valid_ix)));
    outIX = valid_ix(keepIXs);
    
end

if isempty(outIX)
    xx=[];
end

end % of function