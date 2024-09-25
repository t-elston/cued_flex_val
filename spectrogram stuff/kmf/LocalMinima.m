% mins = LocalMinima(x)
%
% finds positions of all strict local minima in input array

function mins = LocalMinima(x)

nPoints = length(x);

Middle = x(2:(nPoints-1));
Left = x(1:(nPoints-2));
Right = x(3:nPoints);

mins = 1+find(Middle < Left & Middle < Right);

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu