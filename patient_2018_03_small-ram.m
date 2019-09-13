%% Patient: 2018_03_Ipsilateral
% Constant Voltage: 3V
% Pulse Width: 90us 
% 1st Stim Block: 140Hz 
% 2nd Stim Block: 20Hz
% 3rd Stim Block: 250Hz
% 4th Stim Block: 70Hz  
clear all; clc; close all;

patient_id = 'Patient-2018-03-Ipsilateral';
stim_names = {'140Hz','20Hz','250Hz','70Hz','__'};

% Patient Path 
patient_path = 'C:\Users\Dakota\dbs-study\patient-2018-03';

% Generate paths to git repository & data 
path_code_files = genpath(patient_path);
path_global_functions = genpath('C:\Users\Dakota\dbs-study\global-functions');
path_data_files = genpath('D:\DBS Study\DBS DATA');

% Add folders & subfolders to path 
addpath(path_code_files);
addpath(path_global_functions);
addpath(path_data_files);

% Load Data and Header Files 
load('D:\DBS Study\DBS DATA\2018_03_Ipsilateral.mat', 'C1');
load('D:\DBS Study\DBS DATA\2018_03_Ipsilateral.mat', 'Header');

Fs = Header.info.sevFs; % sampling frequency
tRaw  = (1:length(C1)) / Fs; % create time vector 


%% Gather Annotations and Times 
commentString = Header.data.Mark.data; % Sting of comments 
commentTimeRaw = Header.data.Mark.time; % time markers per character of comment string
commentMat = strsplit(commentString,"#")'; % splits long string into smaller strings using '#' as a separater 

% Stim ON\OFF Times
startTime = Header.info.startTimeRaw; % start of recording 
commentTime = commentTimeRaw()-startTime; % adjusted time array of comments 

baselineIndex = regexpi(commentString,'\<baseline\>');
baselineTimes = commentTime(baselineIndex);

baselineIndex = regexpi(commentString,'\<baselin\>'); % search for miss-spelled baseline
b_miss_spelled = commentTime(baselineIndex);
baselineTimes = [baselineTimes(1:2) b_miss_spelled baselineTimes(3:end)]; % add miss-spelled baseline to baseline array 


%% Manually specify Begin & End of each block
begin_recording = round(baselineTimes(1)*Fs); % manually specify to aviod artifact
end_recording = round(baselineTimes(5)*Fs); % manually specify to aviod artifact

% capture only necessary data and convert voltage from V to uV 
sig = C1(begin_recording:end_recording)*1e6; 
time = tRaw(begin_recording:end_recording);

clear C1; % Remove unecessary variable 
clear tRaw; % Remove unecessary variable 
%% Identify stim blocks 

% Manually assign begin & end of stim blocks in seconds 
block(1).markers = [baselineTimes(1) baselineTimes(2)];
block(2).markers = [baselineTimes(2) baselineTimes(3)];
block(3).markers = [baselineTimes(3) baselineTimes(4)]; 
block(4).markers = [baselineTimes(4) baselineTimes(5)];

numBlocks = length(block);

% Generate corresponding index values in (samples) and ...
% adjust for time offset when defining indicies 
for i = 1:numBlocks
    block(i).indicies = round((block(i).markers-time(1))*Fs);
    % Indicies cannot be zeros, if this is the case, change to a one 
    if block(i).indicies(1) == 0
        block(1).indicies(1) = 1;
    end
end

%% Add sig & time to blocks structure 
for i = 1: numBlocks
    t1 = block(i).indicies(1);
    t2 = block(i).indicies(2);
    block(i).sig = sig(t1:t2);
    block(i).time = time(t1:t2);
end

clear sig; % Remove unecessary variable 
%% High-Pass Filter
for i = 1:numBlocks
    [block(i).spikeData, filter_params_hp] = filter_high_pass(Fs, block(i).sig);
end

%% Find Stimulation Artifacts 
threshold_art = 400; % threshold just above tallest action potential 
for i = 1:numBlocks
    [block(i).pks_art, block(i).locs_art] = findpeaks(block(i).spikeData,block(i).time,'MinPeakHeight',threshold_art); % MatLab 2018b
end

%% Visualize Blocks 
figure;
for i = 1:4
    subplot(2,2,i)
    plot(block(i).time, block(i).spikeData, block(i).locs_art, block(i).pks_art,'o')
    title(stim_names{i})
end

%% Find Start/Stop of Stimulation
for i = 1:numBlocks
    [block(i).start_art, block(i).stop_art] = start_stop_times(block(i).locs_art);
end

%% Plotting On/Off markers of stimulation

% ii = 1; % *** SELECT A BLOCK TO PLOT ***
% figure(1); clf;
% ax1 = subplot(2,1,1);
% plot(block(ii).time,block(ii).spikeData,block(ii).locs_art,block(ii).pks_art,'o')
% thresh_line = refline(0,threshold_art); % plot threshold line 
% thresh_line.Color = 'r';
% 
% xlim([block(ii).time(1),block(ii).time(end)])
% xlabel('Time(s)')
% ylabel('Voltage(uV)')
% title(stim_names{ii})
% sgtitle('Validate Start/Stop Times of Stimulation')
% 
% ax2 = subplot(2,1,2);
% scale = 500;
% for i = 1:length(block(ii).start_art)
%     line([block(ii).start_art(i) block(ii).start_art(i)], [-scale,scale],'Color','green')
%     line([block(ii).stop_art(i) block(ii).stop_art(i)], [-scale,scale],'Color','red')
%     patch([block(ii).start_art(i), block(ii).stop_art(i), block(ii).stop_art(i), block(ii).start_art(i)], [-scale, -scale, scale, scale], 1, 'facecolor', 'r', 'edgecolor', 'none', 'facealpha', .15  );
% end
% linkaxes([ax1,ax2],'x')

%% Remove Stimulation Artifacts 

artifact_width = 2/1000 * Fs; % Number of samples that span an AP
artifact_array_blank = -artifact_width/2:artifact_width/1.5; % Array with offset to maximize coverage of artifact with little overlap of data

for i = 1:numBlocks
    % subtract initial time value to match with spiikeData indicies 
    artifact_location_index = (block(i).locs_art - block(i).time(1)) * Fs;

    % Matrix of locations of Artifacts 
    artifact_location_matrix = artifact_location_index(:) + artifact_array_blank;

    artifact_matrix = block(i).spikeData(round(artifact_location_matrix)); % Matrix of Artifacts 
    artifact_average = mean(artifact_matrix); % Average stimulation artifact waveform 

    art1D = reshape(artifact_location_matrix,[],1); % create 1D array from matrix of artifact locations 
    artZeros = zeros(size(block(i).spikeData)); % create array of zeros the same length as spikeData 
    artZeros(round(art1D)) = block(i).spikeData(round(art1D)); % create array of artifacts by indexing into spikeData
    
    block(i).cleanSpikeData = block(i).spikeData - artZeros; % subtract artifact from spikeData 
end

%% Visualize Results of Artifact Removal 

% B = 3; % *** SELECT A BLOCK TO PLOT ***
% 
% figure(4); clf;
% ax1 = subplot(2,1,1);
% plot(block(B).time, block(B).spikeData)
% xlim([block(B).time(1),block(B).time(end)])
% title('High-Pass Filtered Signal')
% 
% ax2 = subplot(2,1,2);
% plot(block(B).time, block(B).cleanSpikeData)
% xlim([block(B).time(1),block(B).time(end)])
% title('Remove Artifacts After HP Filtering')
% 
% linkaxes([ax1,ax2],'xy')

%% Determine a quantiative threshold for Action Potentials
threshold = 40; % Choose an appropriate threshold (uV)

% Find Peaks 
for B = 1:numBlocks
    APtime = 2/1000; % duration of an action potential = 2ms
    [block(B).pks_ap_raw, block(B).locs_ap_raw] = findpeaks(block(B).cleanSpikeData,block(B).time,'MinPeakHeight',threshold,'MinPeakDistance',APtime);
end

%% Plot Marked APs
B = 3;
figure(5); clf;
plot(block(B).time,block(B).cleanSpikeData, block(B).locs_ap_raw, block(B).pks_ap_raw,'o')

thresh_line = refline(0,threshold); % plot threshold line 
thresh_line.Color = 'r';

total_APs = [];
for i = 1:numBlocks
    total_APs = cat(2,total_APs, block(i).locs_ap_raw);
end


%% Matrix of APs
spike_width = APtime*Fs; % Number of samples that span an AP
spike_array_blank = -spike_width/2:spike_width/2; % Array centered at zero of length spikeWidth

for i = 1:numBlocks
    locs_remove_ends = block(i).locs_ap_raw(block(i).locs_ap_raw>spike_width & block(i).locs_ap_raw<=(length(block(i).cleanSpikeData)-spike_width)); % remove spikes too close to edges
    locs_time_adjusted = locs_remove_ends - block(i).time(1); % subtract initial time value to match with spiikeData indicies 
    loc_matrix = locs_time_adjusted(:)*Fs + round(spike_array_blank); % Matrix of locations of APs 

    block(i).spikeMat = block(i).cleanSpikeData(round(loc_matrix)); % Matrix of APs 
end

%% Sort Action Potentials
% number of clusters per block
block(1).numClusters = 3; % 140Hz
block(2).numClusters = 3; % 20Hz
block(3).numClusters = 3; % 250Hz
block(4).numClusters = 2; % 70Hz

% Make ALL blocks contain only ONE cluster 
% for ii = 1:numBlocks 
%     block(ii).numClusters = 1;
% end
           
for ii = 1:numBlocks 
    % PCA to cluster spikes 
    [coeff,block(ii).score] = pca(block(ii).spikeMat); % PCA
    Z = linkage(block(ii).score(:,1:2),'ward'); % create a tree of hierarchical clusers using 'ward'
    block(ii).c = cluster(Z,'Maxclust',block(ii).numClusters); % assign spikes to a cluster 
end 
% Plot Clustering results 
spikeColor = [1,0,0; 0,1,0; 0,0,1; 1,0,1; 1,1,0; 0,1,1]; % colors for averages of spikes 
s = 0.8; % color soften coefficient 
stdColor = [1,s,s; s,1,s; s,s,1; 1,s,1; 1,1,s; s,1,1]; % colors for std backgrounds

%% Neuron Subplot 
for B = 1:numBlocks
    figure(B+5); clf;
    for i = 1:block(B).numClusters
        spikeAvg(i,:) = mean(block(B).spikeMat(block(B).c==i,1:end)); % mean action potential
        spikeStd(i,:) = std(block(B).spikeMat(block(B).c==i,1:end)); % std of action potential

        y1 = spikeAvg(i,:)+spikeStd(i,:);
        y2 = spikeAvg(i,:)-spikeStd(i,:);
        Y = [y1 fliplr(y2)];
        X = [spike_array_blank, fliplr(spike_array_blank)];
        
        subplot(3,block(B).numClusters,i)
        fill(X,Y,stdColor(i,:));
        hold on;
        plot(spike_array_blank,spikeAvg(i,:),'Color',spikeColor(i,:),'LineWidth',2);
        xlabel('Samples');ylabel('Voltage (uV)')
        title(['Neuron' num2str(i)])
        ylim([-90,270])
        
    end
    xlim([min(spike_array_blank),max(spike_array_blank)])
    %ylim([min(Y),max(Y)]);
    sgtitle(['Avg AP: ', stim_names{B}])
    
    % Plot marked APs for each cluster 
    subplot(3,block(B).numClusters,block(B).numClusters+1:block(B).numClusters*3)
    plot(block(B).time,block(B).cleanSpikeData)
    hold on
    scatter(block(B).locs_ap_raw, block(B).pks_ap_raw, 10, spikeColor(block(B).c(:),:))
    
end

%% PCA Plot
figure(13); clf;
for B = 1:numBlocks
    subplot(2,2,B);
    scatter3(block(B).score(:,1),block(B).score(:,2),block(B).score(:,3),1,spikeColor(block(B).c(:),:)); % plot each spike based on the first 3 principle components 
    title(stim_names{B});
    xlabel('Component 1');ylabel('Component 2');zlabel('Component 3');
end

sgtitle('PCA Plot')

%% Choose Best AP for each Block 

% Use all APs 
for B = 1:numBlocks 
    block(B).locs_ap = block(B).locs_ap_raw;
    block(B).pks_ap = block(B).pks_ap_raw;
end

% *** Uncomment to select specific neurons ***

block(1).neurons = [3]; % Block 1 = 140Hz
block(2).neurons = [1,3];   % Block 2 = 20Hz
block(3).neurons = [1,3];   % Block 3 = 250Hz
block(4).neurons = [2];     % Block 4 = 70Hz

ylimits = [-200, 500];
figure(15);
for B = 1:numBlocks
    block(B).locs_ap = [];
    block(B).pks_ap = [];
    for i = 1:length(block(B).neurons)
        block(B).locs_ap = cat(2,block(B).locs_ap, block(B).locs_ap_raw(block(B).c==block(B).neurons(i)));
        block(B).pks_ap = cat(2,block(B).pks_ap, block(B).pks_ap_raw(block(B).c==block(B).neurons(i))); 
    end
    subplot(1,numBlocks,B)
    plot(block(B).time,block(B).cleanSpikeData)
    hold on
    scatter(block(B).locs_ap, block(B).pks_ap,10)
    ylim(ylimits)
    title(stim_names{B});
end

%% Firing Rate Plotting 
% for ii = 1:numBlocks
%     AP_locations = block(ii).locs_ap'; % time locations of action potentials 
%     % Create structure of AP locations separated by cluster 
%     for i = 1:block(ii).numClusters
%         myfield = strcat('Neuron',num2str(i));
%         APstruc.(myfield) = AP_locations(block(ii).c==i,1);
%     end
% 
%     seconds_per_bin = 0.1; % length of bins in seconds 
%     bin_num = round(length(time)/(Fs*seconds_per_bin)); % calculated number of bins 
%     scale = 30;
%     figure(7+ii);clf;
%     names = fieldnames(APstruc); % field names in structure
%     for clust_i = 1:block(ii).numClusters
%         subplot(block(ii).numClusters,1,clust_i)
%         h = histogram(APstruc.(names{clust_i}),bin_num); % binned AP location
% 
%         scale = max(h.Values);
%         for a = 1:length(block(ii).start_art) % plot stim patches
%             patch([block(ii).start_art(a), block(ii).stop_art(a), block(ii).stop_art(a), block(ii).start_art(a)], [0, 0, scale, scale], 1, 'facecolor', 'r', 'edgecolor', 'none', 'facealpha', .15  );
%         end
%         title((names{clust_i}))
%         xlabel('Time(s)')
%         ylabel('Firing Rate')
%         ylim([0,max(h.Values)+10])
% 
% 
%     end
%     sgtitle(['Firing Rate: ' stim_names{ii}])
%     
% end



%% LFP: Band Pass Filtering 
for B = 1:numBlocks
    [block(B).lfp, filter_params_bp] = filter_band_pass(Fs, block(B).sig);
end

%% Separate ON/OFF
% *** Define LFP Snippet Size ***
lfp_snip_size = 500/1000; % used to ensure lfp snippet indicies do not go out of bounds
lfp_snip_size_time_offset = lfp_snip_size + time(1);

for B = 1:numBlocks
    lfp_snip_size_time_offset = lfp_snip_size + block(B).time(1);
    block(B).locs_baseline = [];
    block(B).locs_on = [];
    block(B).locs_off = [];
    for i = 1:length(block(B).start_art)
        % BASELINE 
        block(B).locs_baseline = block(B).locs_ap(... 
                                                block(B).locs_ap > lfp_snip_size_time_offset & ...
                                                block(B).locs_ap < block(B).start_art(1));
        % ON Stim
        block(B).locs_on = cat(2, block(B).locs_on, ...
                                           block(B).locs_ap(...
                                           block(B).locs_ap > block(B).start_art(i) & ...
                                           block(B).locs_ap < block(B).stop_art(i))); 
        % OFF Stim
        if i < length(block(B).start_art) 
            block(B).locs_off = cat(2, block(B).locs_off, ...
                                                block(B).locs_ap(...
                                                block(B).locs_ap > block(B).stop_art(i) & ...
                                                block(B).locs_ap < block(B).start_art(i+1))); 
        else % prevent out of bounds indexing (i+1)
            block(B).locs_off = cat(2, block(B).locs_off, ...
                                                block(B).locs_ap(...
                                                block(B).locs_ap > block(B).stop_art(i) & ...
                                                block(B).locs_ap < block(B).time(end)) - lfp_snip_size);
        end
    end
end

%% Validate ON/OFF Separation 
% figure(10);clf;
% for B = 1:numBlocks
%     plot(block(B).locs_baseline,ones(1,length(block(B).locs_baseline)),'b*')
%     hold on
%     plot(block(B).locs_on,ones(1,length(block(B).locs_on)),'r*')
%     plot(block(B).locs_off,ones(1,length(block(B).locs_off)),'g*')
% end

%% LFP Snippets 
for B = 1:numBlocks
    block(B).locs_baseline_index = round((block(B).locs_baseline - block(B).time(1))*Fs);
    block(B).locs_on_index = round((block(B).locs_on - block(B).time(1))*Fs);
    block(B).locs_off_index = round((block(B).locs_off - block(B).time(1))*Fs);
end

lfp_snip_size_samples = lfp_snip_size*Fs; % convert from time (ms) to samples
lfp_array_blank = -lfp_snip_size_samples/2:lfp_snip_size_samples/2; % Array centered at zero of length spikeWidth

%% Create a Matrix of all LFP Snippets 
% Downsampling LFP Data
n_times_upper_cutoff = 10; % sampling frequency = 10 * upper band-pass frequency 
downsampling_factor = round(Fs/(filter_params_bp.Fp2 * n_times_upper_cutoff));

tic
for B = 1:numBlocks
    loc_matrix_lfp = block(B).locs_on_index(:) + round(lfp_array_blank); % Matrix of locations of APs 
    block(B).lfpMat_on = block(B).lfp(round(loc_matrix_lfp)); % Matrix of APs 
    block(B).lfpMat_on_ds = downsample(block(B).lfpMat_on',downsampling_factor)'; % implement downsampling 
    block(B).lfpMat_on = []; % clear non-downsampled data after it is used to preserve memory 
    
    loc_matrix_lfp = block(B).locs_off_index(:) + round(lfp_array_blank); % Matrix of locations of APs 
    block(B).lfpMat_off = block(B).lfp(round(loc_matrix_lfp)); % Matrix of APs 
    block(B).lfpMat_off_ds = downsample(block(B).lfpMat_off',downsampling_factor)';
    block(B).lfpMat_off = [];
    
    loc_matrix_lfp = block(B).locs_baseline_index(:) + round(lfp_array_blank); % Matrix of locations of APs 
    block(B).lfpMat_baseline = block(B).lfp(round(loc_matrix_lfp)); % Matrix of APs 
    block(B).lfpMat_baseline_ds = downsample(block(B).lfpMat_baseline',downsampling_factor)';
    block(B).lfpMat_baseline = [];
end
toc
clear loc_matrix_lfp

%% Create a Matrix of all LFP Snippets 
tic
rows = 3; % Baseline, Off, On
cols = numBlocks; 
y_limits = [-15,15];
colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.4660, 0.6740, 0.1880];

figure(11);clf;
for B = 1:numBlocks
    % Average AP Aligned LFP
    block(B).lfp_avg_on = mean(block(B).lfpMat_on_ds);
    block(B).lfp_avg_off = mean(block(B).lfpMat_off_ds);
    block(B).lfp_avg_baseline = mean(block(B).lfpMat_baseline_ds);
    
    time_array = round(downsample(lfp_array_blank/Fs*1000,downsampling_factor));

    subplot(rows,cols,B)
    plot(time_array, block(B).lfp_avg_baseline,'Color',colors(1,:),'LineWidth',2)
    ylim(y_limits)
    title([stim_names{B}])
    if B==1
        ylabel('Baseline')
    end
    
    subplot(rows,cols,cols+B)
    plot(time_array, block(B).lfp_avg_off,'Color',colors(2,:),'LineWidth',2)
    ylim(y_limits)
    if B==1
        ylabel('Off Stim')
    end
    
    subplot(rows,cols,cols*2+B)
    plot(time_array, block(B).lfp_avg_on,'Color',colors(3,:),'LineWidth',2)
    ylim(y_limits)
    if B==1
        ylabel('On Stim')
    end
end 
sgtitle(['AP Aligned Average LFP: ', patient_id])

%% Plot AP Count Dependence
% AP Count
for B = 1:numBlocks
    block(B).ap_count_on = length(block(B).locs_on);
    block(B).ap_count_off = length(block(B).locs_off);
    block(B).ap_count_baseline = length(block(B).locs_baseline);
end

% RMS
for B = 1:numBlocks
    block(B).rms_on = rms(block(B).lfp_avg_on);
    block(B).rms_off = rms(block(B).lfp_avg_off);
    block(B).rms_baseline = rms(block(B).lfp_avg_baseline);
end

% Create arrays of ap_count and rms for easier plotting 
ap_count = []; 
lfp_rms = [];

for B = 1:numBlocks
    ap_count = cat(2,ap_count,block(B).ap_count_on);
    lfp_rms = cat(2,lfp_rms,block(B).rms_on);
    
    ap_count = cat(2,ap_count,block(B).ap_count_off);
    lfp_rms = cat(2,lfp_rms,block(B).rms_off);
    
    ap_count = cat(2,ap_count,block(B).ap_count_baseline);
    lfp_rms = cat(2,lfp_rms,block(B).rms_baseline);
end

figure(13);clf
scatter(ap_count, lfp_rms,'o')
title('AP Count vs RMS')
xlabel('AP Count')
ylabel('LFP RMS')

%% Power Spectrum of Each LFP Average 
downsampled_Fs = Fs/downsampling_factor;

ylim_low = 0;
ylim_high = 2;

figure(14); clf;

for B = 1:numBlocks
    for on_off_baseline = [1 2 3]
        
        if on_off_baseline == 1
            signal = block(B).lfp_avg_baseline;
            subplot(3,4,B);
            color = colors(1,:);

        elseif on_off_baseline == 2
            signal = block(B).lfp_avg_off;
            subplot(3,4,B + 4);
            color = colors(2,:);
            
        elseif on_off_baseline == 3
            signal = block(B).lfp_avg_on;
            subplot(3,4,B + 8);
            color = colors(3,:);
        end
        
        
        L = length(signal); % length of signal
        NFFT = 2^nextpow2(L)+2; % increase padding 
        f = downsampled_Fs/2*linspace(0,1,NFFT/2+1); % frequency basis for single sided spectrum

        % FFT power spectrum.
        Y = fft(signal,NFFT)/L; % calculate fourier transform
        pxx = abs(Y(1:NFFT/2+1)); % convert Y to signal sided power spectrum
        loglog(f,pxx,'k');
        semilogx(f,abs(Y(1:NFFT/2+1)),'Color',color,'LineWidth',2);
        ylim([ylim_low,ylim_high])
        
        % Titles & Labels 
        if on_off_baseline==1
            title([stim_names{B}])
            ylabel('Baseline')
        elseif on_off_baseline==2 && B==1
            ylabel('Off Stim')
        elseif on_off_baseline==3 && B==1
            ylabel('On Stim')
        end
    end
end

sgtitle(['Power Spectrum of AP aligned LFP: ', patient_id])
 
%% LFP Bootstrapping 

smallest_count = 700; % smallest AP count in a group ***(for bootstrapping)***
epochs = 1000; % number of times to randomly sample 

tic
for B = 1:numBlocks 
    B
    for i = 1:epochs
        rand_rows = randperm(length(block(B).lfpMat_on_ds),smallest_count); % randomly select lfp snippets 
        block(B).lfp_bs_mat_on = block(B).lfpMat_on_ds(rand_rows,:); % create matrix of randomly selected snippets        
        block(B).lfp_bs_avg_on = mean(block(B).lfp_bs_mat_on); 
        
        rand_rows = randperm(length(block(B).lfpMat_off_ds),smallest_count); % randomly select lfp snippets 
        block(B).lfp_bs_mat_off = block(B).lfpMat_off_ds(rand_rows,:); % create matrix of randomly selected snippets        
        block(B).lfp_bs_avg_off = mean(block(B).lfp_bs_mat_off); 
        
        rand_rows = randperm(length(block(B).lfpMat_baseline_ds),smallest_count); % randomly select lfp snippets 
        block(B).lfp_bs_mat_baseline = block(B).lfpMat_baseline_ds(rand_rows,:); % create matrix of randomly selected snippets        
        block(B).lfp_bs_avg_baseline = mean(block(B).lfp_bs_mat_baseline); 
        
    end
end
toc
%% Plotting Bootstrapping 
rows = 3; % Baseline, Off, On

cols = numBlocks; 
y_limits = [-15,15];
colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.4660, 0.6740, 0.1880];

figure(15);clf;
for B = 1:numBlocks    
    time_array = round(downsample(lfp_array_blank/Fs*1000,downsampling_factor));

    subplot(rows,cols,B)
    plot(time_array, block(B).lfp_bs_avg_baseline,'Color',colors(1,:),'LineWidth',2)
    ylim(y_limits)
    title([stim_names{B}])
    if B==1
        ylabel('Baseline')
    end
    
    subplot(rows,cols,cols+B)
    plot(time_array, block(B).lfp_bs_avg_off,'Color',colors(2,:),'LineWidth',2)
    ylim(y_limits)
    if B==1
        ylabel('Off Stim')
    end
    
    subplot(rows,cols,cols*2+B)
    plot(time_array, block(B).lfp_bs_avg_on,'Color',colors(3,:),'LineWidth',2)
    ylim(y_limits)
    if B==1
        ylabel('On Stim')
    end
end 
sgtitle(['AP Aligned LFP: ', patient_id])
%%
a = zeros(10,1);
parfor i = 1:10
    a(i) = i
end
a
%% LFP RMS Bootstrapping
Neuron_array = [1]; % neurons to plot; NOTE: uncomment subplot code 'Neuron-1'
block_array = [2,4,1,3]; % in ascending order of frequency 
colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.4660, 0.6740, 0.1880];
block_name = {'140Hz','20Hz','250Hz','70Hz'};

smallest_count = 700; % smallest AP count in a group ***(for bootstrapping)***
epochs = 1000; % number of times to randomly sample 

countMat = zeros([3,length(blocks)]); % initialize matrix of AP counts
rmsMat = zeros([3,length(blocks)]); % initialize matrix of rms 

rmsAvgmat = zeros(length(block_array),3);
rmsStdmat = zeros(length(block_array),3);

rmscell = cell(length(block_array),3);
powerMat = cell(length(block_array),3);
lfpPower = zeros(1,epochs);
for block = block_array
    for on_off_baseline = [1 2 3]
        lfpMat = lfpcell{block,on_off_baseline};
        
        lfpAvgrms = zeros(1,epochs);
        for i = 1:epochs

            rand_rows = randperm(length(lfpMat(:,1)),smallest_count); % randomly select lfp snippets 
            lfpMat_bootstrap = lfpMat(rand_rows,:); % create matrix of randomly selected snippets 
            
            lfpAvg = mean(lfpMat_bootstrap); 
            lfpAvgrms(i) = std(lfpAvg);
            
            lfpPower(i) = mean(abs(lfpAvg).^2);
            
            
        end
        
        rmscell{block,on_off_baseline} = lfpAvgrms; % matrix of rms at each epoch
        powerMat{block,on_off_baseline} = lfpPower; % matrix of power at each epoch
        
        
        %lfpAvg = mean(lfpMat);
%         rmsAvgmat(block,on_off_baseline) = std(lfpAvg);
% 
%         lfprms = std(lfpMat);
%         rmscell{block,on_off_baseline} = lfprms;
% 
%         rmsStdmat(block,on_off_baseline) = std(rmscell{block,on_off_baseline});
% 

        Block = block
        On_Off_Basleine = on_off_baseline
    end
end
    
%% LFP RMS Plots 
% Initilize matricies 
rmsMat = zeros(epochs,3,length(blocks));
rmsMatmean = zeros(length(blocks),3);
errMat = zeros(length(blocks),3);

rmsPC = zeros(epochs,3,length(blocks));
rmsPCmean = zeros(length(blocks),3); % percent change of rms from baseline
errMatPC = zeros(length(blocks),3);
fontsize = 14;
for block = block_array
    
    rms_baseline = mean(rmscell{block,3}(:)); % baseline rms for each epoch 
   
    for on_off_baseline = [1 2 3]
        errMat(block,on_off_baseline) = std(rmscell{block,on_off_baseline});
        rmsMatmean(block,on_off_baseline) = mean(rmscell{block,on_off_baseline});
        rmsMat(:,on_off_baseline,block) = rmscell{block,on_off_baseline}(:);
        
        
%       rms_percent_change = (rmscell{block,on_off_baseline}'-rms_baseline)./rms_baseline*100; 
        baseline_norm = rmscell{block,on_off_baseline}'-rms_baseline;
        rms_percent_change = baseline_norm*100; 
        
        rmsPC(:,on_off_baseline,block) = rms_percent_change;
        rmsPCmean(block,on_off_baseline)  = mean(rms_percent_change);
        errMatPC(block,on_off_baseline)  = std(rms_percent_change);
    end
end

xx = [1 2 3];
figure(41);clf;
% Plot RMS Values 
i = 1;
for block = block_array
        
    subplot(1,length(block_array),i)
    i = i+1;
    errorbar(xx,fliplr(rmsMatmean(block,:)),fliplr(errMat(block,:)),'o-')
    %plot(xx,fliplr(rmsMat(block,:)),'o-')
    hold on
    title([block_name{block} ' Stimulation'],'fontsize',fontsize/1.5)
    xlim([0.5,3.5])
    ylim([1.5,8])
    ylabel('Voltage (uV)','fontsize',fontsize/1.5)
    xticks([1 2 3])
    xticklabels({'Baseline','OFF','ON'})
        
end
sgtitle(['Ipsilateral P#1: LFP Average RMS using Bootstrapping, AP Count: ', num2str(smallest_count),', Epochs: ', num2str(epochs)],'fontsize',fontsize)

xx = [1 2 3];
figure(42);clf;
% Percent Change from baseline 
i = 1;
for block = block_array
        
    subplot(1,length(block_array),i)
    i = i+1;
    errorbar(xx,fliplr(rmsPCmean(block,:)),fliplr(errMatPC(block,:)),'o-')
    %plot(xx,fliplr(rmsMat(block,:)),'o-')
    hold on
    title([block_name{block} ' Stimulation'],'fontsize',fontsize/1.5)
    xlim([0.5,3.5])
    ylim([-150,300])
    ylabel('Percent Cahnge (%)','fontsize',fontsize/1.5)
    xticks([1 2 3])
    xticklabels({'Baseline','OFF','ON'})
        
end
sgtitle(['Ipsilateral P#1: Percent Change of LFP Average RMS using Bootstrapping, AP Count: ', num2str(smallest_count),', Epochs: ', num2str(epochs)],'fontsize',fontsize)

% Bar Chart of Percent Change from Baseline
figure(43);clf;
%rmsPCmean = rmsMatmean;
%errMatPC = errMat;
%errMatPC = errMatPC;
i = 1;
for block = block_array
    
    subplot(1,length(block_array),i)
    i = i+1;
    
    bar(xx,fliplr(rmsPCmean(block,:)))
    hold on
    er = errorbar(xx,fliplr(rmsPCmean(block,:)),fliplr(errMatPC(block,:)));
    er.LineStyle = 'none';  
    er.Color = [0 0 0];
    title([block_name{block} ' Stimulation'],'fontsize',fontsize/1.5)
    xlim([0.5,3.5])
    ylim([-200,300])
    %ylim([0,12])
    ylabel('Percent Change (%)','fontsize',fontsize/1.5)
    xticks([1 2 3])
    xticklabels({'Baseline','OFF','ON'})
    
        
end
sgtitle(['Ipsilateral P#1: Percent Change of LFP Average RMS using Bootstrapping, AP Count: ', num2str(smallest_count),', Epochs: ', num2str(epochs)],'fontsize',fontsize)




%% LFP POWER
rmsMat = zeros(epochs,3,length(blocks));
rmsMatmean = zeros(length(blocks),3);
errMat = zeros(length(blocks),3);

for block = block_array
    for on_off_baseline = [1 2 3]
        errMat(block,on_off_baseline) = std(powerMat{block,on_off_baseline});
        rmsMatmean(block,on_off_baseline) = mean(powerMat{block,on_off_baseline});
        rmsMat(:,on_off_baseline,block) = powerMat{block,on_off_baseline}(:);
    end
end

xx = [1 2 3];
figure(46);clf;
% Plot RMS Values 
i = 1;
for block = block_array
        
        subplot(1,length(block_array),i)
        i = i+1;
        errorbar(xx,fliplr(rmsMatmean(block,:)),fliplr(errMat(block,:)),'o-')
        %plot(xx,fliplr(rmsMat(block,:)),'o-')
        hold on
        title([block_name{block} ' Stimulation'],'fontsize',fontsize/1.5)
        xlim([0.5,3.5])
        ylim([0,75])
        ylabel('Mean AP Aligned LFP RMS Voltage (uV)','fontsize',fontsize/1.5)
        xticks([1 2 3])
        xticklabels({'Baseline','OFF','ON'})
        
end
sgtitle(['Ipsilateral P#1: LFP Average RMS using Bootstrapping, AP Count: ', num2str(smallest_count),', Epochs: ', num2str(epochs)],'fontsize',fontsize)

%% Check Normal Distribution 
xlimit1 = -400;
xlimit2 = 400;
num_bin = 30;
figure(50);clf;
for block = block_array
    for on_off_baseline = [1 2 3]
        %data2hist = rmsPC(:,on_off_baseline,block) - mean(rmsPC(:,on_off_baseline,block));
        data2hist = rmsPC(:,on_off_baseline,block);
        
        if block == 1
            subplot(length(block_array),3,on_off_baseline)
            histogram(data2hist,num_bin)
            xlabel('Percent Change (%)')
            ylabel('Count')
            %xlim([xlimit1,xlimit2])
            if on_off_baseline ==1
                title('Block:1, ON')
            elseif on_off_baseline ==2
                title('Block:1, OFF')
            elseif on_off_baseline ==3
                title('Block:1, Baseline')
            end
        elseif block == 2
            subplot(length(block_array),3,on_off_baseline+3)
            histogram(data2hist,num_bin)
            xlabel('Percent Change (%)')
            ylabel('Count')
            %xlim([xlimit1,xlimit2])
            if on_off_baseline ==1
                title('Block:2, ON')
            elseif on_off_baseline ==2
                title('Block:2, OFF')
            elseif on_off_baseline ==3
                title('Block:2, Baseline')
            end
        elseif block == 3
            subplot(length(block_array),3,on_off_baseline+6)
            histogram(data2hist,num_bin)
            xlabel('Percent Change (%)')
            ylabel('Count')
            %xlim([xlimit1,xlimit2])
            if on_off_baseline ==1
                title('Block:3, ON')
            elseif on_off_baseline ==2
                title('Block:3, OFF')
            elseif on_off_baseline ==3
                title('Block:3, Baseline')
            end
        elseif block == 4
            subplot(length(block_array),3,on_off_baseline+9)
            histogram(data2hist,num_bin)
            xlabel('Percent Change (%)')
            ylabel('Count')
            %xlim([xlimit1,xlimit2])
            if on_off_baseline ==1
                title('Block:4, ON')
            elseif on_off_baseline ==2
                title('Block:4, OFF')
            elseif on_off_baseline ==3
                title('Block:4, Baseline')
            end
        end
    end
end
sgtitle('Ipsilateral P#1: Histograms of RMS Percent Change from Baseline')
%% Statistical Analysis 
block_name = {'140Hz','20Hz','250Hz','70Hz'};

rmsPCstats = zeros(3,2,4);
for block = block_array
    [p,tbl,stats] = anova1(rmsPC(:,:,block),{'ON','OFF','Baseline'},'off');
    figure(51+block);clf;
    [c,m] = multcompare(stats,'Alpha',0.01);
    block
    c
    rmsPCstats(:,:,block) = m;
    title([block_name(block)])
end

%% RMS PC Stats plot

figure(55);clf;
i = 1;
fontsize = 14;
for block = block_array
    
    rms_data = flipud(rmsPCstats(:,1,block)); % flip rms array to plot baseline first, then OFF, then ON
    err_data = flipud(rmsPCstats(:,2,block)); % array for error bars: Standard Error (multiply by 'sqrt(epochs)' for SD 
    
    subplot(1,length(block_array),i)
    i = i+1;
    bar(xx,rms_data)
    hold on
    er = errorbar(xx,rms_data,err_data*2); % multiply by 2 to easily visualize significance 
    er.LineStyle = 'none'; % remove connecting lines from error bar chart
    er.Color = [0 0 0]; % color of error bars  
    
    title([block_name{block} ' Stimulation'],'fontsize',fontsize/1.4)
    xlim([0.5,3.5])
    ylim([-60,120])
    %ylim([0,12])
    ylabel('Percent Change (%)','fontsize',fontsize/1.4)
    xticks([1 2 3])
    xticklabels({'Baseline','OFF','ON'})
    
        
end
sgtitle(['Ipsilateral P#1: Percent Change of LFP Average RMS using Bootstrapping, AP Count: ', num2str(smallest_count),', Epochs: ', num2str(epochs)],'fontsize',fontsize/1.2)





