% Out = NormVector(In)
%
% takes a matrix that consists of column vectors and normalized each one.

function Out = NormVector(In)

nDim = size(In, 1);

Norms = sqrt(sum(abs(In).^2, 1));

Out = In ./ repmat(Norms,nDim,1);

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu