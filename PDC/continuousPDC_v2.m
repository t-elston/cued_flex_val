% routine which produces windowed PDC estimates, similar to a coherence
% routine

function [runningCh1toCh2,runningCh2toCh1,permuted_runningCh1toCh2,permuted_runningCh2toCh1] = continuousPDC_v2(ch1,ch2,window_size,overlap,decimationFactor,info_crit)



% inputs: ch1,ch2,window_size, overlap, decimation factor, information
% criteria
% outputs: time-resolved PDC estimates for 
% ch1 -->ch2 
% ch2--> ch1,
% permuted ch1-->ch2
% permuted ch2 -->ch1




runningCh1toCh2=[];
runningCh2toCh1=[];
permuted_runningCh1toCh2=[];
permuted_runningCh2toCh1=[];
freq_band=[4:12];

% step 1

% decimate data
overlap = overlap/decimationFactor;
window_size = window_size/decimationFactor;
LFP1(:,1) = decimate(ch1,decimationFactor);
LFP2(:,1) = decimate(ch2,decimationFactor);

% step 2
% load data into single array
data = horzcat(LFP1,LFP2);
clear LFP1;
clear LFP2;
clear ch1;
clear ch2;


% step 3
% start crawling through the data
len = length(data(:,1));

for i = 1:window_size:len
    
    % step 4
    % find the part of the data you want to do the PDC analysis on
    window_start = i-(.5*overlap);
    window_end = i+(.5*overlap);
    if window_end > len
    else
    currentPDCblock=[];
    currentPDCblock=data(i:i+overlap,:);
    
    % step 5
    % generate permutated data of the same time series
    shuffle_len = length(currentPDCblock(:,1));

     % create random permutation of that order
        rand_idx = randperm(shuffle_len);

      % shuffle the LFP data based on the randomized index
      shuffled_ch1 = currentPDCblock(rand_idx,1);
      current_shuffledPDC1_block = horzcat(shuffled_ch1,currentPDCblock(:,2));
      shuffled_ch2 = currentPDCblock(rand_idx,2);
      current_shuffledPDC2_block = horzcat(currentPDCblock(:,1),shuffled_ch2);
      
      
    % step 6
    % carry out the PDC on the respective blocks
      blockPDC=[];
      shuffledPDC1=[];
      shuffledPDC2=[];
      blockPDC     = PDC(currentPDCblock, info_crit, freq_band);
      shuffledPDC1 = PDC(current_shuffledPDC1_block, info_crit, freq_band);
      shuffledPDC2 = PDC(current_shuffledPDC2_block, info_crit, freq_band);
  
      current_ch1_to_ch2 = sum(nanmean(blockPDC(:,:,2)));
      current_ch2_to_ch1 = sum(nanmean(blockPDC(:,:,1)));
      current_shuffled_ch1_to_ch2 = sum(nanmean(shuffledPDC1(:,:,2)));
      current_shuffled_ch2_to_ch1 = sum(nanmean(shuffledPDC2(:,:,1)));
      
      % step 7
      % add PDCS to output arrays
      runningCh1toCh2            = vertcat(runningCh1toCh2,current_ch1_to_ch2);
      runningCh2toCh1            = vertcat(runningCh2toCh1,current_ch2_to_ch1);
      permuted_runningCh1toCh2   = vertcat(permuted_runningCh1toCh2,current_shuffled_ch1_to_ch2);
      permuted_runningCh2toCh1   = vertcat(permuted_runningCh2toCh1,current_shuffled_ch2_to_ch1);
      
      % step 8
      % clear variables for next cycle
      shuffled_ch1=[];
      shuffled_ch2=[];
      current_shuffledPDC1_block=[];
      current_shuffledPDC2_block=[];
      rand_idx=[];
      blockPDC=[];
      current_ch1_to_ch2=[];
      current_ch2_to_ch1=[];
      current_shuffled_ch1_to_ch2=[];
      current_shuffled_ch2_to_ch1=[];
      
    end % of if loop
end % of for loop

      return
      
      
      
  

