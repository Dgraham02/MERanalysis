% Global Function
% Purpose: Remove stimualtion artifacts from data 
% Notes:

function [cleanSpikeData, artAvg] = remove_stim_artifacts(general_data, post_hp_data, locs_art_remove)

% Define variables from structures
Fs = general_data.Fs;
time = general_data.time;
spikeData = post_hp_data.spikeData;

locs_art1 = locs_art_remove; % select removal spikes 
artTime = 1/1000; % length of artifact
artWidth = artTime*Fs; % Number of samples that span an AP
artLen = -artWidth/2:artWidth/1.8; % Array with offset to maximize coverage of artifact with little overlap of data
locs_art = locs_art1 - time(1); % subtract initial time value to match with spiikeData indicies 
loc_artMat = locs_art(:)*Fs + artLen; % Matrix of locations of Artifacts 

artMat = spikeData(round(loc_artMat)); % Matrix of Artifacts 
artAvg = mean(artMat); % Average stimulation artifact waveform 

art1D = reshape(loc_artMat,[],1); % create 1D array from matrix of artifact locations 
artZeros = zeros(size(spikeData)); % create array of zeros the same length as spikeData 
artZeros(round(art1D)) = spikeData(round(art1D)); % create array of artifacts by indexing into spikeData
cleanSpikeData = spikeData - artZeros; % subtract artifact from spikeData 
end