% CrossTabPlot(X1, X2) - from raw values
% or CrossTabPlot(m) - from cross-table
%
% Draws a plot showing actual and predicted values under the
% independence hypothesis in a cross-table.
%
% See also CrossTab

function CrossTabPlot(X1, X2)


if (nargin>1)
	[Table, Vals1, Vals2] = CrossTab(X1, X2);
else
	Table = X1;
	Vals1 = 1:size(Table, 1);
	Vals2 = 1:size(Table, 2);
end

n = sum(Table(:));

p1 = sum(Table, 1) / n;
p2 = sum(Table, 2) / n;

Predicted = n*(p2 * p1);

v1 = length(Vals1);
v2 = length(Vals2);

clf
for x1 = 1:v1
	for x2=1:v2
		subplot(v1+1, v2+1, 1+x2 + x1*(v2+1))
		tit = sprintf('%s,%s', num2str(Vals1(x1)), num2str(Vals2(x2)));
		hold on
		h1 = bar(1, Table(x1,x2), 'b');
		h2 = bar(2, Predicted(x1,x2), 'r');
		hold off
		set(gca, 'XTick', []);
		xlim([.5 2.5])
		% legend
		if (x1==1 & x2==v2)
			legend('Actual', 'Predicted');
		end
	end
end

ForAllSubplots(sprintf('ylim([0 %d])', ceil(max([Table(:) ; Predicted(:)]))));

% marginals
for x1=1:v1
	subplot(v1+1, v2+1, 1+x1*(v2+1));
	bar(sum(Table(x1,:),2));
	ylim([0 max(sum(Table,2))]);
%	set(gca, 'XTickLabel', 'Total');
	ylabel(['X1 = ', num2str(Vals1(x1))]);
end

for x2=1:v2
	subplot(v1+1, v2+1, 1+x2);
	bar(sum(Table(:,x2),1));
	ylim([0 max(sum(Table,1))]);
%	set(gca, 'XTickLabel', 'Total');
	title(['X2 = ', num2str(Vals2(x2))]);
end




% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu