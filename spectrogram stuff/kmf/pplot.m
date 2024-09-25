function pplot(X)
% pplot(X) - plots multiple graphs with a pause between each one

MaxX = size(X, 1);
MinX = 1;

MaxY = max(X(:));
MinY = min(X(:));

axis([MinX MaxX MinY MaxY]);

for i=1:size(X,2)
	h = plot(X(:,i));
	%axis([MinX MaxX MinY MaxY]);
	pause;
	delete(h);
end;

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu