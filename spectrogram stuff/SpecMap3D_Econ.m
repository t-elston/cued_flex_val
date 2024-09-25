 function [mapCohT mapCh1T mapCh2T mapCohLG mapCh1LG mapCh2LG mapCohMG mapCh1MG mapCh2MG mapCohHG mapCh1HG mapCh2HG newX newY newXYts XI YI XCohT XCh1T XCh2T XCohLG XCh1LG XCh2LG XCohMG XCh1MG XCh2MG XCohHG XCh1HG XCh2HG yg fg] = SpecMap3D_Econ(ch1,ch2,downsampling,nFFT,fs,wlength,overlap,tapers,detrend,keep,Series,PositionX,PositionY,PositionY_ts,EventsCommentStopping_Recording);
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
thetaband=[4 12];
lowgamma=[35 55];
midgamma=[60 80];
highgamma=[90 110];

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
ch1PSD=yo(:,:,1,1);
ch2PSD=yo(:,:,2,2);

%set frequency band of interest to analyse
the=nFFT/fs*(thetaband(1,1));
ta=nFFT/fs*(thetaband(1,2));

lo=nFFT/fs*(lowgamma(1,1));
wgamma=nFFT/fs*(lowgamma(1,2));

mi=nFFT/fs*(midgamma(1,1));
dgamma=nFFT/fs*(midgamma(1,2));

hig=nFFT/fs*(highgamma(1,1));
hgamma=nFFT/fs*(highgamma(1,2));

%perform averaging across frequency band of interest
thetaCOH=coh(the:ta,:);
mean_thetaCOH=mean(thetaCOH);

thetaCh1=ch1PSD(the:ta,:);
mean_thetaCh1=mean(thetaCh1);

thetaCh2=ch2PSD(the:ta,:);
mean_thetaCh2=mean(thetaCh2);

lgammaCOH=coh(lo:wgamma,:);
mean_lgammaCOH=mean(lgammaCOH);

lgammaCh1=ch1PSD(lo:wgamma,:);
mean_lgammaCh1=mean(lgammaCh1);

lgammaCh2=ch2PSD(lo:wgamma,:);
mean_lgammaCh2=mean(lgammaCh2);

mgammaCOH=coh(mi:dgamma,:);
mean_mgammaCOH=mean(mgammaCOH);

mgammaCh1=ch1PSD(mi:dgamma,:);
mean_mgammaCh1=mean(mgammaCh1);

mgammaCh2=ch2PSD(mi:dgamma,:);
mean_mgammaCh2=mean(mgammaCh2);

hgammaCOH=coh(hig:hgamma,:);
mean_hgammaCOH=mean(hgammaCOH);

hgammaCh1=ch1PSD(hig:hgamma,:);
mean_hgammaCh1=mean(hgammaCh1);

hgammaCh2=ch2PSD(hig:hgamma,:);
mean_hgammaCh2=mean(hgammaCh2);

%Correct tracking; not really necessary anymore since correction is USUALLY
%done beforehand in the larger script.

%[sPX sPY]=Fix_tracking(Series,PositionX, PositionY, PositionY_ts,0);
sPX=PositionX;
sPY=PositionY;
sts=PositionY_ts;

%This part cuts out the OF session at the end if there is one. Or,
%rather, this only processes the first part of the recording session
finish=find(round(PositionY_ts)==round(EventsCommentStopping_Recording(1,1)),1);

sPX=sPX(1:finish);
sPY=sPY(1:finish);
sts=PositionY_ts(1:finish);

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
nPositionY_ts=(0:length(newXYts)-1)'/25;
mapCoh=zeros(length(nPositionY_ts),1);

for i=1:length(nPositionY_ts);
    dist=(to-nPositionY_ts(i)).^2;
    [m,ind]=min(dist);
        mapCohT(i)=mean_thetaCOH(ind(1));
    mapCh1T(i)=mean_thetaCh1(ind(1));
    mapCh2T(i)=mean_thetaCh2(ind(1));
        mapCohLG(i)=mean_lgammaCOH(ind(1));
    mapCh1LG(i)=mean_lgammaCh1(ind(1));
    mapCh2LG(i)=mean_lgammaCh2(ind(1));
        mapCohMG(i)=mean_mgammaCOH(ind(1));
    mapCh1MG(i)=mean_mgammaCh1(ind(1));
    mapCh2MG(i)=mean_mgammaCh2(ind(1));
        mapCohHG(i)=mean_hgammaCOH(ind(1));
    mapCh1HG(i)=mean_hgammaCh1(ind(1));
    mapCh2HG(i)=mean_hgammaCh1(ind(1));
end

%Plot spectra at each downsampled time point
% figure;
% plot3c(newX,newY,newXYts,mapCoh);

%binning and averaging the 3-D scatter data into matrix form for plotting
%purposes
[XI,YI]=meshgrid(-300:10:300,0:10:420);
ZCohT = bin2mat(newX,newY,mapCohT,XI,YI); 
[XCohT] = nanmoving_average2(ZCohT,1,1);
ZCh1T = bin2mat(newX,newY,mapCh1T,XI,YI);
[XCh1T] = nanmoving_average2(ZCh1T,1,1);
ZCh2T = bin2mat(newX,newY,mapCh2T,XI,YI);
[XCh2T] = nanmoving_average2(ZCh2T,1,1);

ZCohLG = bin2mat(newX,newY,mapCohLG,XI,YI); 
[XCohLG] = nanmoving_average2(ZCohLG,1,1);
ZCh1LG = bin2mat(newX,newY,mapCh1LG,XI,YI);
[XCh1LG] = nanmoving_average2(ZCh1LG,1,1);
ZCh2LG = bin2mat(newX,newY,mapCh2LG,XI,YI);
[XCh2LG] = nanmoving_average2(ZCh2LG,1,1);

ZCohMG = bin2mat(newX,newY,mapCohMG,XI,YI); 
[XCohMG] = nanmoving_average2(ZCohMG,1,1);
ZCh1MG = bin2mat(newX,newY,mapCh1MG,XI,YI);
[XCh1MG] = nanmoving_average2(ZCh1MG,1,1);
ZCh2MG = bin2mat(newX,newY,mapCh2MG,XI,YI);
[XCh2MG] = nanmoving_average2(ZCh2MG,1,1);

ZCohHG = bin2mat(newX,newY,mapCohHG,XI,YI); 
[XCohHG] = nanmoving_average2(ZCohHG,1,1);
ZCh1HG = bin2mat(newX,newY,mapCh1HG,XI,YI);
[XCh1HG] = nanmoving_average2(ZCh1HG,1,1);
ZCh2HG = bin2mat(newX,newY,mapCh2HG,XI,YI);
[XCh2HG] = nanmoving_average2(ZCh2HG,1,1);

%Plot the scaled/averaged spectral map
% figure;
% imagesc(XCoh);
% axis xy;


end

