 %% Snippets of LFP

% on = 1, off = 2, baseline = 3
function [lfpMat,lfpAvg,spikeAvg] = LFP_Snip_Func(lfp_data,time,LFPstruc,Fs,c,locs2,spikeMat,block,on_off_baseline,lfp_snippet_size,Neuron)
    ii = block;
    i = on_off_baseline;
    
    adjust_window = 60; % adjust for being too close to the start point

    
    on_block_field = strcat('on_block',num2str(ii));
    off_block_field = strcat('off_block',num2str(ii));
    baseline_block_field = strcat('baseline_block',num2str(ii));

    if i==1
        lfp_locs = LFPstruc.(on_block_field);
    elseif i==2
        lfp_locs = LFPstruc.(off_block_field);
    elseif i==3
        lfp_locs = LFPstruc.(baseline_block_field);
    end

    lfpTime = lfp_snippet_size/1000; % seconds of lfp
    lfpWidth = lfpTime*Fs; % Number of samples that span an AP
    lfpLen = -lfpWidth/2:lfpWidth/2; % Array centered at zero of length spikeWidth

    c_begin = find(locs2==lfp_locs(2)); % index of first cluster value 
    c_end = find(locs2==lfp_locs(end)); % index of last cluster value 

    locs3 = locs2(c_begin:c_end)*Fs; % locations of spikes within specified region
    locMat1 = locs3(:) + round(lfpLen); % matrix of locations 

    % 
    locMat2 = locMat1(find(locMat1(:,ceil(end/2))>(lfp_locs*Fs)-20 & locMat1(:,ceil(end/2))<(lfp_locs*Fs)+20));
    locMat = (locMat2 + round(lfpLen))-time(1)*Fs;

    locMat_mid = locMat1(:,ceil(end/2));
    lfp_c = sum(locMat_mid>((lfp_locs*Fs)-20) & locMat_mid<((lfp_locs*Fs)+20),2);

    CC = c(cast(lfp_c,'logical')); % turn cluster array back to logical 
    if ii==1
        CC = CC(adjust_window:end); % adjust cluster number 
        locMat = locMat(adjust_window:end,:); % adjust spike number 
    end

    lfpMat = lfp_data(round(locMat(CC==Neuron,1:end))); % Matrix of LFP for cluster 'n'

    %lfpSamples = linspace(-lfpWidth/2,lfpWidth/2,length(lfpLen)); % samples points of lfp snippet to orient at zero
    %spikeSamples = linspace(-spikeWidth/2,spikeWidth/2,length(spikeLen));  % samples points of AP snippet to orient at zero

    % Plotting Parameters
    %lfp_snip_time = lfpSamples/Fs*1000; % used for labeling X-ticks on plot
    %spike_snip_time = spikeSamples/Fs*1000; % used for labeling X-ticks on plot

    spikeAvg = mean(spikeMat(CC==Neuron,1:end)); % average action potential waveform
    lfpAvg = mean(lfpMat); % mean lfp waveform
    % lfpAvg = lfpMat(50,:); % ***VALIDATION***
    lfpStd = std(lfpMat); % std of lfp waveform
end