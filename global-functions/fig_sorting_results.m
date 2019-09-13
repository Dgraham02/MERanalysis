% Global Function
% Purpose: Plot results of spike sorting 
% Notes: 
function [f1,f2,f3] = fig_sorting_results(general_data, post_hp_data,clean_data, ap_data,save_fig)

patient_id = general_data.patient_id;
Fs = general_data.Fs;
time = general_data.time;
start_art = post_hp_data.start_art;
stop_art = post_hp_data.stop_art;
numBlocks = post_hp_data.numBlocks;
cleanSpikeData = clean_data.cleanSpikeData;
spikeMat = ap_data.spikeMat;
pks = ap_data.pks;
locs = ap_data.locs;
spikeLen = ap_data.spikeLen;
numClusters = ap_data.numClusters;
c = ap_data.c;


% Colors for plotting
spikeColor = [1,0,0; 0,1,0; 0,0,1; 1,0,1; 1,1,0; 0,1,1]; % colors for averages of spikes 
s = 0.8; % color soften coefficient 
stdColor = [1,s,s; s,1,s; s,s,1; 1,s,1; 1,1,s; s,1,1]; % colors for std backgrounds

% Plot peak locations of different clusters 
f1 = figure('Name','Peak Locations');clf;
plot(time,cleanSpikeData)
hold on
for i = 1:numClusters
    plot(locs(c==i)+time(1),pks(c==i),'o','Color',spikeColor(i,:))
end
legend('Filtered/Clean Signal', 'Neuron 1', 'Neuron 2','Location','northwest')
title('Verifying PCA Clustering: Peak Locations')
xlabel('Time(s)')
ylabel('Voltage(uV)')

if save_fig == true
    fig_path = strcat(patient_path,'\figures\',f.Name,'_',patient_id);
    savefig(f,fig_path)
end

%% Neuron Subplot 
f2 = figure('Name','AP Snippets'); clf;

for i = 1:numClusters

    spikeAvg(i,:) = mean(spikeMat(c==i,1:end)); % mean action potential
    spikeStd(i,:) = std(spikeMat(c==i,1:end)); % std of action potential
    
    y1 = spikeAvg(i,:)+spikeStd(i,:);
    y2 = spikeAvg(i,:)-spikeStd(i,:);
    Y = [y1 fliplr(y2)];
    X = [spikeLen, fliplr(spikeLen)];
    subplot(1,numClusters,i)
    fill(X,Y,stdColor(i,:));
    hold on;
    plot(spikeLen,spikeAvg(i,:),'Color',spikeColor(i,:),'LineWidth',2);
    xlabel('Samples');ylabel('Voltage (uV)')
    title(['Neuron' num2str(i)])
    ylim([-90,270])
    
end
xlim([min(spikeLen),max(spikeLen)])
%ylim([min(Y),max(Y)]);
sgtitle('Average AP Snippet: Ipsilateral Patient #1')

if save_fig == true
    fig_path = strcat(patient_path,'\figures\',f.Name,'_',patient_id);
    savefig(f,fig_path)
end


%% Firing Rate Plotting 
AP_locations = locs'; % time locations of action potentials 
% Create structure of AP locations separated by cluster 
for i = 1:numClusters
    myfield = strcat('Neuron',num2str(i));
    APstruc.(myfield) = AP_locations(c==i,1);
end

freq_names = {'140Hz','20Hz','250Hz','70Hz','50Hz'};

seconds_per_bin = 0.1; % length of bins in seconds 
bin_num = round(length(time)/(Fs*seconds_per_bin)); % calculated number of bins 

f3 = figure('Name','Firing Rate');clf;
names = fieldnames(APstruc); % field names in structure
for ii = 1:numClusters
    subplot(numClusters,1,ii)
    h = histogram(APstruc.(names{ii}),bin_num); % binned AP location
    
    scale = max(h.Values);
    for i = 1:length(start_art) % plot stim patches
        patch([start_art(i), stop_art(i), stop_art(i), start_art(i)], [0, 0, scale, scale], 1, 'facecolor', 'r', 'edgecolor', 'none', 'facealpha', .15  );
    end
    title((names{ii}))
    xlabel('Time(s)')
    ylabel('Firing Rate')
    ylim([0,max(h.Values)+10])
    
    for f_name = 1:numBlocks
        dimstep = (f_name-1)/6.6;
        dim = [.24+dimstep .55 .1 .3];
        str = freq_names(f_name);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
    end
    for f_name = 1:numBlocks
        dimstep = (f_name-1)/6.6;
        dim = [.24+dimstep .1 .1 .3];
        str = freq_names(f_name);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
    end
    
end
sgtitle('Firing Rate of Neurons: Ipsilateral Patient #1')

if save_fig == true
    fig_path = strcat(patient_path,'\figures\',f.Name,'_',patient_id);
    savefig(f,fig_path)
end

end