% test PDC


%v2 by TW Elston on 20/Feb/2017
% adds AIC model order calculation
% doubles amount of noise in signal.


%v3 by TW Elston on 27/Feb/2017
% adds and plots permutations to noise
% generates frequency bands of interest in a more discreet manner


% 
randmix=false;  % make true to randomly mix items in waveforms
randnoise=true; % make true to add noise to waveforms
sn=10;           %signal to noise for randnoise

Fs=480;         % sampling freq of data
freqHz1=8;      % freq of generated test sine
freqHz2=8;      % freq of generated test sine
upbound=12;     % freq band limits in Hz
lobound=1;
offset1=0;      % phase shift for wave 1 in samples
offset2=-10;     % phase delay wave2  (-ve values) in samples
Fc1=10;          % create 10 cycles of test waveform
Fc2=10;          % create 10 cycles of test waveform

% check which wave is slower, scale wave parameters to that
if freqHz1 < freqHz2
Fc2 = Fc1*(freqHz2/freqHz1);
elseif freqHz1 > freqHz2     
Fc1 = Fc2*(freqHz1/freqHz2);
end



momax=20;       % maximum order of fit for PDC;
model_order = [];

period1=Fs/freqHz1;       % find period of sine
period2=Fs/freqHz2;       

freq_band = lobound:upbound;

bandstart=Fs/lobound;   % set boundarys of test in samples
bandend=Fs/upbound;

% maybe the comb shape occurs because it's checking so regularly
freq = Fs./freq_band;


hz=Fs./freq-1;          % convert from samples to Hz for display



wave1=makesine(period1,offset1,Fc1);   % make test waves
wave2=makesine(period2,offset2,Fc2);
%wave1=wave1+1;   % don't do this
%wave2=wave2+1;

if randmix
  wave1=wave1(randperm(length(wave1)));
  wave2=wave2(randperm(length(wave2)));
end

wave1_block=[];
wave2_block=[];
noise_wave1=[];
noise_wave2=[];
num_perms=100;

if randnoise
    for perm_num = 1:num_perms 
  noise_wave1 = awgn(wave1,sn,'measured');
         pause_length = rand;
         pause(pause_length);
  noise_wave2 = awgn(wave2,sn,'measured');
  
  wave1_block(perm_num,:) = noise_wave1;
  wave2_block(perm_num,:) = noise_wave2;
  
  
    end % of for
end % of if 
model_orderWaves=[];
model_orderWaves(1,:,1) = wave1_block(1,:);
model_orderWaves(2,:,1) = wave2_block(1,:);

% calculate univresal model order
[~,~,moAIC,~] = tsdata_to_infocrit(model_orderWaves,momax);

 model_order = moAIC;

   waveboth=[];
   ch2_leads=[];
   ch1_leads=[];
   testwave1=[];
   testwave2=[];
   
% now do PDC stuff
for perm_num = 1:num_perms 
    waveboth=[];
    testwave1=[];
    testwave2=[];
    testwave1 = wave1_block(perm_num,:).';
    testwave2 = wave2_block(perm_num,:).';
    
    waveboth= cat(2,testwave1,testwave2);   % cat two waves into one array
    pdc=[];
    [pdc] = PDC(waveboth,model_order,freq);
    val1= pdc(:,1,1);       % ch i  to itself
    val2= pdc(:,1,2);       % ch i to j
    val3= pdc(:,2,1);
    val4= pdc(:,2,2);
    
    ch1_leads(:,perm_num) = val2;
    ch2_leads(:,perm_num) = val3;
end
    
ch1_leads_mean = mean(ch1_leads,2);
ch2_leads_mean = mean(ch2_leads,2);

all_perms = cat(2,ch1_leads,ch2_leads);

 maxval=1;
 minval=0;

figure;
subplot(2,1,1)
h1 = plot(hz,ch1_leads(:,1),'b');
ylim([minval maxval])
hold on
h2=plot(hz,ch1_leads,'b');
ylim([minval maxval])
h3=plot(hz,ch1_leads_mean,'r','lineWidth',2);
ylim([minval maxval])
hold off
legend([h1 h3],{'Permutations','Mean'})
title('Ch1 Leads')

subplot(2,1,2)
ylim([minval maxval])
h1 = plot(hz,ch2_leads(:,1),'b');
hold on
h2=plot(hz,ch2_leads,'b');
ylim([minval maxval])
h3=plot(hz,ch2_leads_mean,'r','lineWidth',2);
ylim([minval maxval])
hold off
legend([h1 h3],{'Permutations','Mean'})
title('Ch2 Leads')
