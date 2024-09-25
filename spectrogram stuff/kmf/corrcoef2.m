function xy = corrcoef2(x,y)
%CORRCOEF2 Correlation coefficients for 2 multivariate data sets
%   CORRCOEF2(X,Y), where X is n by m and Y is n by p gives a
%   m by p correlation matrix out
%   
%   If C is the covariance matrix, C = COV(X), then CORRCOEF(X) is
%   the matrix whose (i,j)'th element is
%
%          C(i,j)/SQRT(C(i,i)*C(j,j)).
%
%   See also COV2, CORRCOEF

c = cov2(x,y);
vx = var(x);
vy = var(y);

xy = c./sqrt(vx'*vy);


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu