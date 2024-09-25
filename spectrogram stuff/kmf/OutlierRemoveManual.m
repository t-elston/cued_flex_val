function [Included, Excluded] = OutlierRemove(Data)
% [Included Excluded] = OutlierRemove(Data)
%
% An interactive, iterative function to remove
% outliers from multivariate data, based on the chi^2
% plot

% Parameter: number worst to display.

nWorst = 8;

nIncluded = size(Data,1);
nDims = size(Data,2);

Included = (1:nIncluded);
Excluded = [];

while(1)
	% Calculate Mahalanobis distances
	d = Data(Included, 1:nDims);
	m = mahal(d,d);
	
	%Do Chi-Squared plot
	chi2plot(m, nDims);

	[Sorted, Index] = sort(m);
	
	% calculate covariance matrix
	CovMat = cov(d);
	fprintf('Covariance trace %g det %g\n', trace(CovMat), det(CovMat));
	
	fprintf('%d Worst outliers:\n', nWorst);
	for i = Index(nIncluded:-1:nIncluded+1-nWorst)'
		fprintf('point %d with Mahalanobis distance %f, p = %f\n',...
				Included(i), m(i), 1-chi2cdf(m(i), nDims));
	end;
	
	s = input('Input number of top values to remove, or n to stop.> ', 's');
	
	if (s == 'n')
		break;	
	else
		nToRemove = str2num(s);
		ToExclude = Included(Index(nIncluded + 1 - nToRemove : nIncluded));
		Excluded = [Excluded(:)' , ToExclude(:)'];
		Included = setdiff(Included, ToExclude);
		nIncluded = nIncluded - nToRemove;
	end;
end;


	
	

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu