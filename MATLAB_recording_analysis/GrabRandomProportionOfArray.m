function [outIX,leftOverIX] = GrabRandomProportionOfArray(indata,proportion)

% this function works on one 1-D vector at a time
% Thomas Elston 
% telston@nurhopsi.org
% May 5 2021

outIX=[];

randIX = randperm(numel(indata));
keepIXs = randIX(1:floor(proportion*numel(indata)));
outIX = indata(keepIXs);

leftOverIX = randIX(floor(proportion*numel(indata)):end);
LeftOverOutIX = indata(leftOverIX);


end % of function