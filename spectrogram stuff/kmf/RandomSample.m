% v = RandomSample(m, n)
%
% produces a vector of m integers in the range 1:n
% without replacement.

function v = RandomSample(m, n)

rp = randperm(n);

v = rp(1:m);


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu