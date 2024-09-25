%%Written by Erin Boykin, 2011
% v02     by Thomas Elston, 2017
%————————————————————————————————-
%Inputs: X – time series for analysis (rows are time points, columns are
% channels)
% pmax – max order of AR model fit (perform an AIC or BIC analysis prior
% in order to make sure max is set sufficiently high)
% freq – vector of discrete frequencies at which gc is calculated
% doGPDC - boolean of whether to do GPDC
%
%Outputs: pdc – pdc values in form frequency x channel i x channel j, where
% causality is from channel i to channel j
%————————————————————————–
%The ARfit package for computing AR models is available here : http://www.gps.caltech.edu/~tapio/arfit/
% it is required in order for this function to run

function [pdc] = PDC_v02(X,pmax,freq,doGPDC)

Nc = size(X,2); %% Nc = channels

pdc = zeros(length(freq),Nc,Nc);

[w,B,Sig]=arfit(X,1,pmax,'fpe','zero');

% Thom added this
if isnan (Sig)
    pdc = NaN;
    return
end

Bf = getBf(Nc,B,freq);

for xx = 1:Nc
for yy = 1:Nc
if (yy ~= xx)
for ff = 1:length(freq)
pdc(ff,yy,xx) = abs(Bf(xx,yy,ff))/sqrt(Bf(:,yy,ff)'*Bf(:,yy,ff));

% GPDC 
if doGPDC
 pdc(ff,yy,xx) = (abs(Bf(xx,yy,ff))*(1/sqrt(Sig(xx,xx))))/abs(sqrt(Bfnorm'*Bf(:,yy,ff)));
end % of whether to do GPDC

end
end
end
end

function [Bf] = getBf(m,B,freq)

popt = size(B,2) / m;
A0 = eye(m);
Bf = zeros(m,m,length(freq));

for g=1:length(freq)
for ch=1:m
Bj = B(:,[ch:m:m*popt]);
for j = 1:popt
Bf(:,ch,g) = Bf(:,ch,g) + Bj(:,j).*exp(1i*j*freq(g));
end
end
Bf(:,:,g) = A0 - Bf(:,:,g);
end

