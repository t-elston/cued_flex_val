function a = clip(b, min, max)
% clip(b, min, max)
% takes a matrix and replaces any elements below min
% with min and any above max with max

shape = size(b);
b = b(:);
b(find(b<min)) = min;
b(find(b>max)) = max;
a = reshape(b, shape);



% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu