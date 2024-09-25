% [a b R2 U V] = CanonCor(X, Y)
%
% does a canonical correlation analysis on two sets of data X and Y
%
% R2 is the returned vector of cannonical covariates
% a and b are the projections that produce them
%
% Each column of a or b corresponds to a projection
% with the first column being the projection onto the
% 1st cannonical variable, etc.
%
% U and V are the cannonical variables themselves.
% Each row corresponds to a case number and each
% column to a cannonical variable.

function [a, b, R2, U, V] = CanonCor(X, Y)

% Make covariance matrices
XSize = size(X, 2);
YSize = size(Y, 2);

BigCov = cov([X, Y]);
CXX = BigCov(1:XSize, 1:XSize);
CYY = BigCov(XSize+1:end, XSize+1:end);
CXY = BigCov(1:XSize, XSize+1:end);

% make -1/2 powers ...

CXXMH = CXX ^ -0.5;
CYYMH = CYY ^ -0.5;

% matrix to do svd on ...
M = CXXMH * CXY * CYYMH;

% do svd
[c s d] = svd(M);

% then calculate a and b

a = CXXMH * c;
b = CYYMH *d;
R2 = diag(s);

if (nargout > 3)
	U = X*a;
	V = Y*b;
end;

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu