% Local Function 
% Purpose: Filter raw signal to extract action potentials 
% Notes: Using filter Order=1 to reduce aliasing of stim artifacts 

function [lfp, filter_params_bp] = filter_band_pass(Fs, sig)

% Filter Parameters 
filter_params_bp.Order = 10;
filter_params_bp.Fp1 = 3; 
filter_params_bp.Fp2 = 55;

% Filter Data
x = detrend(double(sig)); % detrend signal 

% Filter Parameters
Order  = filter_params_bp.Order; % Order
Fp1 = filter_params_bp.Fp1; % 1st pass frequency 
Fp2 = filter_params_bp.Fp2; % 2nd pass frequency 

% Create filter 
h = fdesign.bandpass('N,F3dB1,F3dB2',Order,Fp1,Fp2,Fs);
Hd = design(h, 'butter');

% Apply Filter 
lfp_1 = filtfilt(Hd.sosmatrix, Hd.scalevalues, x); 

lfp = lfp_1'; % transpose lfp data to correct dimensions

end