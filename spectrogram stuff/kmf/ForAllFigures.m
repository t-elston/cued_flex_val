% ForAllFigures(Command)
%
% Selects every figure
% and executes Command.

function ForAllFigures(Command);

for a=get(0, 'Children')';
	set(0, 'CurrentFigure', a);
	eval(Command);
end;

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu