% Local Function 
% Purpose: Filter raw signal to extract action potentials 
% Notes: Using filter Order=1 to reduce aliasing of stim artifacts 

function [spikeData, filter_params_hp] = high_pass_filter(Fs, sig)

% Filter Parameters 
filter_params_hp.Order = 1;
filter_params_hp.Fc = 300; 

x = detrend(double(sig)); % detrend signal 

% Filter Parameters
N  = filter_params_hp.Order; % Order
Fc = filter_params_hp.Fc; % Cutoff Frequency

% Create filter 
h  = fdesign.highpass('N,F3dB', N, Fc, Fs); 
Hd = design(h, 'butter');

% Apply Filter 
spikeData = filtfilt(Hd.sosmatrix, Hd.scalevalues, x); 

end