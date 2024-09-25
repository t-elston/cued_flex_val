% [Table, Vals1, Vals2] = CrossTab(X1, X2)
%
% A slightly different and hopefully faster way of 
% crosstabulating than statistics toolbox's crosstab.
%
% X and Y can be vectors of anything - don't need to be
% integers 1...n.
% 
% Table is a matrix giving counts that a particular 
% value of X and Y occur together.
%
% Vals1 and Vals2 give the actual value for the 
% nth row or column of Table.

function [Table, Vals1, Vals2] = CrossTab(X1, X2)

if (length(X1) ~= length(X2))
	error('X1 and X2 should be the same length');
end;

Vals1 = unique(X1);
nVals1 = length(Vals1);

Vals2 = unique(X2); 
nVals2 = length(Vals2);

Table = zeros(nVals1, nVals2);

for i1 = 1:nVals1
	for i2 = 1:nVals2
		Table(i1, i2) = sum(X1 == Vals1(i1) & X2 == Vals2(i2));
	end;
end;

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu