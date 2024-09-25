function [Included, Excluded] = OutlierRemove(Data, pVal)
% [Included Excluded] = OutlierRemove(Data, pVal)
%
% An automatic, iterative function to remove
% outliers from multivariate data, At each step, the
% worst outlier will be removed if the p-value for
% its Mahalanobis distance is less than pVal.

nIncluded = size(Data,1);
nDims = size(Data,2);

Included = (1:nIncluded);
Excluded = [];

while(1)
	% Calculate Mahalanobis distances
	d = Data(Included, 1:nDims);
	m = mahal(d,d);
	
	%Do Chi-Squared plot
	%chi2plot(m, nDims);
	%pause;

	[Worst, Index] = max(m);
	ToExclude = Included(Index);
	p = 1-chi2cdf(Worst, nDims);
	
	fprintf('Worst outlier: ');
	fprintf('point %d with Mahalanobis distance %f, p = %f\n', ToExclude, Worst, p);

	if (p > pVal)
		break;	
	else
		Excluded = [Excluded(:)' , ToExclude];
		Included = setdiff(Included, ToExclude);
		nIncluded = nIncluded - 1;
	end;
end;

% final Chi^2 plot
% chi2plot(m, nDims);

	
	

% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu