% Global Function 
% Purpose: Create a matrix of AP snippets 
% Notes: 
function [spikeMat, locs, spikeLen] = make_ap_matrix(clean_data, general_data)

Fs = general_data.Fs;
time = general_data.time;
cleanSpikeData = clean_data.cleanSpikeData;
APtime = clean_data.APtime;
locs1 = clean_data.locs1;


spikeWidth = APtime*Fs; % Number of samples that span an AP
spikeLen = -spikeWidth/2:spikeWidth/2; % Array centered at zero of length spikeWidth

locs2 = locs1(locs1>spikeWidth & locs1<=(length(cleanSpikeData)-spikeWidth)); % remove spikes too close to edges
locs = locs2 - time(1); % subtract initial time value to match with spiikeData indicies 
locMat = locs(:)*Fs + round(spikeLen); % Matrix of locations of APs 

spikeMat = cleanSpikeData(round(locMat)); % Matrix of APs 

%whitenoise = rms(cleanSpikeData)*randn(49000,1); % random noise 
%noiseMat = reshape(whitenoise,[1000,49]);
%spikeMat = [spikeMat;noiseMat];
end