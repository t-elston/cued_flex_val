% PlotLine(w, b, Linespec)
%
% Plots a line w'*x + b = 0
% it's a big line so use 'axis manual' before running
% this function
%
% LineSpec is optional

function PlotLine(w,b)



% if w = [0 0] cant do it
if (w(1) == 0 & w(2) == 0)
	return;
end;

% see if slope is close to horizontal or close to vertical
if(abs(w(2))>abs(w(1)))
	% close to horizontal
	x0 = -1e4;
	x1 = 1e4;
	y0 = -(b+w(1)*x0)/w(2);
	y1 =-(b+w(1)*x1)/w(2);
else
	% close to vertical
	y0 = -1e4;
	y1 = 1e4;
	x0 = -(b+w(2)*y0)/w(1);
	x1 = -(b+w(2)*y1)/w(1);
end;

if (nargin<3) plot([x0 x1], [y0 y1]);
else plot([x0 x1], [y0 y1], LineSpec);
end;

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu