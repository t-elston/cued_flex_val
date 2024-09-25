% Quantile(X, y)
%
% Computes the quantiles that show how the values in y are distributed 
% relative to the distribution of X.  So if X was a sample from N(0,1),
% and y was 0, the answer would be 0.5.

function q = Quantile(X, y)

% remove Nans from X
X(find(isnan(X))) = [];

nX = length(X);
ny = length(y);

ToSort = [X(:) ; y(:)];

Label = [ones(nX,1) ; zeros(ny,1)];

[Sorted Index] = sort(ToSort);

CumDist = cumsum(Label(Index));

qAll(Index) = CumDist/nX;

q = qAll(nX+1:end);

% don't allow NaNs through
q(find(isnan(y))) = NaN;

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu