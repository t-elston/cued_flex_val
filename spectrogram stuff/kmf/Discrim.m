% [Disc y] = Discrim(x, g)
%
% Performs a multiple linear discriminant analysis on data x that belong
% to groups g.  x is a matrix whose rows correspond to cases
% and whose columns correspond to variables.  g is an array with
% one entry per case containing a positive integer marking the group
% for the corresponding case.
%
% Disc gives the discriminant coefficents; y gives the points in discriminant space

function [Disc, y] = Discrim(x, g);

nGroups = max(g);
[nCases nDim] = size(x);

Means = zeros(nGroups, nDim);

for i=1:nGroups
	Means(i,:) = mean(x(find(g==i),:),1);
end

Between = cov(Means);

Within = cov(x - Means(g,:));

[v d] = eig(Within \ Between);

[Sorted Order] = sort(diag(d));

Order2 = flipud(Order);

Disc = v(:,Order2);

y = x * Disc;
