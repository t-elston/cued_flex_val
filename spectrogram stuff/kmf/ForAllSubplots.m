% ForAllSubplots(Command)
%
% Selects every subplot of the current figure
% and executes Command.

function ForAllSubplots(Command);

for a=get(gcf, 'Children')';
	set(gcf, 'CurrentAxes', a);
	eval(Command);
end;

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu