function [notch_filts] = make_notch_filts(fs)

% make notch filters: 3 @ 60Hz

Nfilts = 3; % number filters: 60, 120, 180
notch_filts = cell(1,Nfilts);

for f = 1:Nfilts
    
    freq = f*60;
    
    d = designfilt('bandstopfir', ... response type
        'FilterOrder',1000, ... filter order
        'CutoffFrequency1',freq-2, ... frequency constraints
        'CutoffFrequency2',freq+2, ...
        'SampleRate',fs); ... sampling rate
        
    notch_filts{f} = d;
    
end

end