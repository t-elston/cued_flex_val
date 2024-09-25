% QuickRegress(x,y)
%
% does a quick scatterplot of y vs. x, plots a lsline
% and prints stats

function QuickRegress(x,y)

nPoints = length(x);

plot(x,y,'*')
lsline

[b bint r rint stats] = regress(y, [x(:), ones(nPoints,1)]);

str = sprintf('y = %fx + %f.  R^2 = %f.  p = %f', b, stats([1 3]));

title(str);

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu