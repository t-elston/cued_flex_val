%%Written by Erin Boykin, 2011
%————————————————————————————————-
%Inputs: X – time series for analysis (rows are time points, columns are
% channels)
% pmax – max order of AR model fit (perform an AIC or BIC analysis prior
% in order to make sure max is set sufficiently high)
% freq – vector of discrete frequencies at which gc is calculated
%
%Outputs: gpdc – gpdc values in form frequency x channel i x channel j, where
% causality is from channel i to channel j
%————————————————————————–
%The ARfit package for computing AR models is available here : http://www.gps.caltech.edu/~tapio/arfit/
% it is required in order for this function to run

function [gpdc] = GPDC(X,pmax,freq)

Nc = size(X,2); %% Nc = channels

gpdc = zeros(length(freq),Nc,Nc);

[w,B,Sig]=arfit(X,1,pmax,'fpe','zero');
Bf = getBf(Nc,B,freq);


for xx = 1:Nc
for yy = 1:Nc
if (yy ~= xx)
for ff = 1:length(freq)
gpdc(ff,yy,xx) = (abs(Bf(xx,yy,ff))*(1/sqrt(Sig(xx,xx))))/abs(sqrt(norm(Bf)'*Bf(:,yy,ff)));
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