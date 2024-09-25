% PlotConic(M, v, c, LineSpec, nPoints)
%
% Plots the conic section that is the solution
% of x'*M*x + v'*x + c = 0
%
% --- under most circumstances (i.e. if M is not singular)
%
% LineSpec defaults to blue
% nPoints defaults to 100
%
% returns handle to plot

function h = PlotConic(M, v, c, LineSpec, nPoints)

if (nargin<4) LineSpec = 'b'; end;
if (nargin<5) nPoints = 100; end;

v = v(:);

Minv = inv(M);

Center = -Minv*v/2;

Const = v'*Minv*v/4 - c;

% if (Const<0) warning('something funny going on'); end;

Th = (0:nPoints)*2*pi/nPoints;

A = M(1,1);
B = M(1,2) + M(2,1);
C = M(2,2);

r = sqrt (Const./(A*cos(Th).^2 + B*cos(Th).*sin(Th) + C*sin(Th).^2));
Minr = min(r);

% if(~isreal(r)) warning('not all r points are real'); end;

x(1,:) = r.*cos(Th) + Center(1,ones(1, nPoints+1));
x(2,:) = r.*sin(Th) + Center(2,ones(nPoints+1,1));

%plot(x(1,:), x(2,:));
h = []; % h is handle vector
hold on
for i=1:nPoints
	if (isreal(x(:,i:i+1)) & isfinite(x(:,i:i+1)))
		h = [h, plot(x(1,i:i+1), x(2,i:i+1), LineSpec)];
		%scatter(x(1,i), x(2,i));
	end;
end;


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu