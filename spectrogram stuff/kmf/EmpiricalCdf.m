% [x c] = EmpiricalCdf(data);
%
% computes the cdf of a set of input data. 
% this is a very simple function that basically only involves a sort

function [x, c] = EmpiricalCdf(data);

nData = length(data);

[x order] = sort(data);

c = (1:nData)'/nData;


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu