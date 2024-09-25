function [coh ch1PSD ch2PSD mean_thetaCOH mean_thetaCh1 mean_thetaCh2 newX newY newXYts XI YI XCoh XCh1 XCh2 yg fg to] = SpecMap3D(ch1,ch2,downsampling,nFFT,fs,wlength,overlap,tapers,detrend,keep,thetaband,sPX,sPY,sts);
%In most contexts, ch1=PFC and ch2= HPC
%This script takes Neuroexplorer data and computes multi-tapered spectra to
%plot onto the figure8 maze
%Usual usage,
%downsampling=4;
%nFFT=1600;
%fs=800;
%wlength=1600;
%overlap=1500;
%tapers=3;
%detrend=1;
%keep=5;
%thetaband=[4 12];
%[mapCoh mapCh1 mapCh2 newX newY newXYts XI YI XCoh XCh1 XCh2 fg yg] = SpecMap3D(Ch2,Ch4,4,1600,800,1600,1500,3,1,5,[4 12],'E',PositionX,PositionY,PositionY_ts,EventsCommentStopping_Recording);
x=ch1;
y=ch2;
%this is the downsampling factor
factor=downsampling;

%downsampling
dx=decimate(x,factor);
dy=decimate(y,factor);

%re-organize data and perform multi-taper spectra calculations
data=[dx,dy];
fs=fs/factor;
nFFT=nFFT/factor;
wlength=wlength/factor;
overlap=overlap/factor;


%figure;
[yg, fg]=mtchd(data,nFFT,fs,wlength,overlap,tapers,detrend,keep);
%PlotMatrix(fg,yg);
[yo,fo,to]=mtchg(data,nFFT,fs,wlength,overlap,tapers,detrend,keep);



%extract spectra of interest from 4-D array
coh=yo(:,:,1,2);
ch1PSD=(yo(:,:,1,1));
ch2PSD=(yo(:,:,2,2));

% % % comment this out when doing many-frequencies
% set frequency band of interest to analyse
the=nFFT/fs*(thetaband(1,1));
ta=nFFT/fs*(thetaband(1,2));

%perform averaging across frequency band of interest
thetaCOH=coh(the:ta,:);
mean_thetaCOH=log10(mean(thetaCOH));

thetaCh1=ch1PSD(the:ta,:);
mean_thetaCh1=log10(mean(thetaCh1));

thetaCh2=ch2PSD(the:ta,:);
mean_thetaCh2=log10(mean(thetaCh2));

%Correct tracking
% [sPX sPY]=Fix_tracking(Series,PositionX, PositionY, PositionY_ts,0);






%downsampling tracking data to match the time points in the spectra
%analyses
newX=decimate(sPX,ceil(length(sPX)/length(to)));
newY=decimate(sPY,ceil(length(sPX)/length(to)));
newXYts=decimate(sts,ceil(length(sPX)/length(to)));

% %median filter and interpolation of the downsampled tracking data
% medY=medfilt1(newY,12);
% medX=medfilt1(newX,12);

%vectorizing timestamps to find what the power and coherence are for the
%frequency band of interest at each X,Y coordinate
nPositionY_ts=(0:length(newXYts)-1)'/50;
mapCoh=zeros(length(nPositionY_ts),1);

for i=1:length(nPositionY_ts);
    dist=(to-nPositionY_ts(i)).^2;
    [m,ind]=min(dist);
    mapCoh(i)=mean_thetaCOH(ind(1));
    mapCh1(i)=mean_thetaCh1(ind(1));
    mapCh2(i)=mean_thetaCh2(ind(1));   
end

%Plot spectra at each downsampled time point
% figure;
% plot3c(newX,newY,newXYts,mapCoh);

%binning and averaging the 3-D scatter data into matrix form for plotting
%purposes
 [XI,YI]=meshgrid(10:5:800,30:5:700);
ZCoh = bin2mat(newX,newY,mapCoh,XI,YI); 
[XCoh] = nanmoving_average2(ZCoh,1,1);
ZCh1 = bin2mat(newX,newY,mapCh1,XI,YI);
[XCh1] = nanmoving_average2(ZCh1,1,1);
ZCh2 = bin2mat(newX,newY,mapCh2,XI,YI);
[XCh2] = nanmoving_average2(ZCh2,1,1);

% XCh1=smooth2a(XCh1,1);
% XCh2=smooth2a(XCh2,1);
% XCoh=smooth2a(XCoh,1);


%Plot the scaled/averaged spectral map
 %figure;
 % imagesc(XCh1);
% view([-90 90]);


end

