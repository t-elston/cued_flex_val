% Clustogram(Fet, Clu, Colors)
%
% Draws a shaded clustogram from the Feature file and cluster file.
% Also Fet must be a 2 column array
% Colors gives the base color of each cluster.  It is optional, and
% will default to [red green blue yellow magenta cyan white]
%

function Clustogram(Fet, Clu, Colors)

if (nargin < 3 | isempty(Colors))
	Colors = [1 0 0 ; 0 1 0 ; 0 0 1 ; 1 1 0 ; 1 0 1; 1 1 0 ; 1 1 1];
end;

Gamma = 2.5;

MinX = min(Fet(:,1));
MaxX = max(Fet(:,1));
MinY = min(Fet(:,2));
MaxY = max(Fet(:,2));

BinsX = MaxX-MinX+1;
BinsY = MaxY-MinY+1;
Im = zeros([BinsY, BinsX, 3]);

% compute bins
BinX = 1+Fet(:,1)-MinX;
BinY = 1+Fet(:,2)-MinY;

% add up histo
for i=1:size(Fet, 1)
	%for c=1:3
		Im(BinY(i), BinX(i),:) = Im(BinY(i), BinX(i),:) + reshape(Colors(Clu(i),:), [1 1 3]);
	%end;
end;

% normalize histo
MaxHist = max(max(max(Im)));

image([MinX, MaxX], [MinY, MaxY], (Im/MaxHist).^(1/Gamma));

% now draw an ellipse
hold on;

Clu2Points = Fet(find(Clu==2), :);

Center=mean(Clu2Points, 1);
Center2 = mean(Fet, 1);
CovMat = cov(Clu2Points);
RootMat = CovMat^0.5;
nPoints = 30;
Th = 0:2*pi/nPoints:2*pi;
Ellipse = 2*RootMat*[cos(Th) ; sin(Th)];
Ellipse = Ellipse + Center(ones(size(Th)), :)';

plot(Ellipse(1,:), Ellipse(2,:), 'w');
hold off;



% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu