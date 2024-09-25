%%Updated by Erin Boykin, 2011
%%Adapted from code by Rosenblum and Pikovsky:
%%http://www.agnld.uni-potsdam.de/~mros/damoco.html
%————————————————————————————————-
%Inputs: y – time series for analysis (rows are time points, columns are
% channels)
% This function is for bivariate series with *2* channels
% Fsamp – sampling frequency of signal
%
%Outputs: c12 – pdm value from time series 1 to time series 2
% c21 – pdm value from time series 2 to time series 1
%————————————————————————–
%This function assumes that theh signals in y are narrow-band enough that
%they have a well-defined phase. Since this project deals with neuronal
%oscillators which are narrow-band, this assumption is upheld.

function [c12 c21] = PDM_2(y,Fsamp)

if size(y,2) ~= 2
error('Must have bivarate time series');
end

npt=length(y);
ncut=1; % number of points to be thrown away in the beginning
% and int the end of signals
npforw = 100; % this is the only parameter for dirc routine
% this parameter corresponds to tau in the paper
%%%% (shift in points to determine delta phi)
% reasonable choice is npforw =

%%%%%%%%%%%%%%%%%%%
%%%%%% Hilbert transform and computation of phases
%%%%%% Phase plots for check
ht=hilbert(y(:,1)); phi1=angle(ht); phi1=phi1(ncut:npt-ncut);
ht=hilbert(y(:,2)); phi2=angle(ht); phi2=phi2(ncut:npt-ncut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi1=unwrap(phi1); phi2=unwrap(phi2);% unwrap phases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Check synchronization index
synchroindex=(mean(cos(phi1-phi2)))^2 + (mean(sin(phi1-phi2)))^2;
if synchroindex > 0.6
disp('WARNING: synchronization index is rather high!!!');
synchroindex
disp('The results might be not significant!');
disp('Comparison of different algorithms is recommended');
end

%%%%%%%% directionality indices
%%%%%%%% EMA %%%%%%%%%%%%%%%%%
[c12 c21]=dirc(phi1,phi2,npforw); % inputs are unwrapped phases

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c12 c21]=dirc(phi1,phi2,npforw) %%% inputs are unwrapped phases

pi2=pi+pi;
npt=length(phi1)-npforw;
dphi1=phi1(npforw+1:npforw+npt);
dphi2=phi2(npforw+1:npforw+npt);
phi1=phi1(1:npt); phi2=phi2(1:npt);
dphi1=dphi1-phi1; dphi2=dphi2-phi2;
phi1=mod(phi1,pi2); phi2=mod(phi2,pi2);

% design matrix is common for both fittings
X = [ones(size(phi1)) sin(phi1) cos(phi1) sin(phi2) cos(phi2)...
sin(phi1-phi2) cos(phi1-phi2) sin(phi1+phi2) cos(phi1+phi2)...
sin(2*phi1) cos(2*phi1) sin(2*phi2) cos(2*phi2)...
sin(3*phi1) cos(3*phi1) sin(3*phi2) cos(3*phi2)];

a1=X\dphi1;
a2=X\dphi2;

c21=sqrt(dzdy2(a1));
c12=sqrt(dzdx2(a2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=dzdx2(a)
r= 4*a(10)^2 + 4*a(11)^2 + 9*a(14)^2 + 9*a(15)^2 + a(2)^2 ...
+ a(3)^2 + a(6)^2 + a(7)^2 + a(8)^2 + a(9)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=dzdy2(a)
r= 4*a(12)^2 + 4*a(13)^2 + 9*a(16)^2 + 9*a(17)^2 + a(4)^2 ...
+ a(5)^2 + a(6)^2 + a(7)^2 + a(8)^2 + a(9)^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%