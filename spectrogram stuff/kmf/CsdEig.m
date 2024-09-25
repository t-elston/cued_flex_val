% CsdEig(Csd, f, nEigs, AnchorChannel, Gamma)
%
% Does a display based on eigen-analysis of a Csd.
% The input Csd should be a 3d array, with the last
% two dimensions being channel number
%
% f is an array of frequencies, as produced by mtcsd.
% default is 0...1
%
% nEigs is the number of eigenvalues to consider
%
% AnchorChannel is the channel which phases are measured
% relative to.
%
% gamma is a contrast modifier for the Complex Color plots
% (see ComplexImage).  Default 0.8

function CsdEig(Csd, f, nEigs, AnchorChannel, Gamma)

if (nargin<2 | isempty(f))	f = ((1:size(Csd,1)) - 1) / (size(Csd,1)); end
if (nargin<3 | isempty(nEigs)) nEigs = 1; end;
if (nargin<4 | isempty(AnchorChannel)) AnchorChannel = 1; end
if (nargin<5 | isempty(Gamma)) Gamma = 0.8; end;

nCh = size(Csd,3);
nFreqs = length(f);
nSubPlots = nEigs+1;


% Eigenvalues of CSD
Ratio = zeros(nFreqs,nEigs);
EigStore = zeros(nFreqs, nCh, nEigs);
	
for f0=1:length(f)
	CsMat = squeeze(Csd(f0,:,:));
	[v d] = eig(CsMat);

%	[MaxEig Index] = max(diag(d));plot(f, real(Ratio2), 'c');

%	Ratio(f0) = MaxEig / sum(diag(d));
%	RotMaxEV = v(:,Index).*abs(v(AnchorChannel,Index))./v(AnchorChannel,Index);
%	EigStore(f0,:) = RotMaxEV';	
	
	% calculate second eigenvalue just for a laugh
	[Sorted Order] = sort(diag(d));
	Index = Order(nCh:-1:(nCh-nEigs+1));
	Ratio(f0,:) = Sorted(nCh:-1:(nCh-nEigs+1))'/sum(diag(d));
	AnchorPhase = v(AnchorChannel,Index) ./ abs(v(AnchorChannel,Index));
	RotMaxEV = v(:,Index) ./ repmat(AnchorPhase,nCh,1);
	EigStore(f0,:,:) = RotMaxEV;	

end

%Plot top eigenvalue ratio
subplot(nSubPlots,1,1)
cla
plot(f, real(Ratio));
axis([0 max(f) 0 1]);
hold on
plot([0 max(f)], [1/nCh,1/nCh], '--');
	
% Plot eigenvectors
for e=1:nEigs
	subplot(nSubPlots,1,e+1);
	ComplexImage(f,1:nCh,EigStore(:,:,e)',Gamma) 
end	


% Written by Kenneth D. Harris 
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu