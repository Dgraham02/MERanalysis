% Global Function 
% Purpose: Band-Pass filter to extract LFP
% Notes: 
function [lfp, Fp1, Fp2] = band_pass_filter(general_data,show_filter)

Fs = general_data.Fs;
sig = general_data.sig; 

Data=sig; % use all data

% Filter Data
x = detrend(double(Data)); % detrend signal 
% Filter Parameters
N  = 10; % Order
Fp1 = 3; % Lower frequency bound
Fp2 = 55; % Upper frequency bound 

% Create filter 
h = fdesign.bandpass('N,F3dB1,F3dB2',N,Fp1,Fp2,Fs);
Hd = design(h, 'butter');

% Apply Filter
lfp1 = filtfilt(Hd.sosmatrix, Hd.scalevalues, x); % implement filter
lfp = lfp1'; % transpose spikeData 

if show_filter == true
fvtool(Hd)
xlim([0,0.1])
ylim([-100,10])
title(['LFP - Band-Passed at ', num2str(Fp1),'-',num2str(Fp2),'Hz'])
end
end