% Global Function
% Putpose: Plot a snippet of data showing complex burst firing of neurons 
% Notes: 
function [f] = fig_complex_burst(general_data, post_hp_data, chosenBaselineTimes,save_fig)

% Define variables from structures 
patient_path = general_data.patient_path;
patient_id = general_data.patient_id;
stim_names = general_data.stim_names;
Fs = general_data.Fs;
time = general_data.time;
spikeData = post_hp_data.spikeData;

ylim_low = min(spikeData); % upper voltage bound 
ylim_high = max(spikeData); % lower voltage bound 

duration = 20000; % length of window (miliseconds)

f = figure('Name','Complex Burst'); clf;

for i = 1:length(chosenBaselineTimes)
    t1 = chosenBaselineTimes(i) - time(1);
    t2 = t1 + duration/1000;
    
    s1 = round(t1*Fs)+1;
    s2 = round(t2*Fs);
    
    sig_window = spikeData(s1:s2); 
    time_window = (1:length(sig_window))/Fs + t1+time(1);
    
    subplot(1,length(chosenBaselineTimes),i)
    plot(time_window,sig_window)
    xlim([time_window(1),time_window(end)])
    xlabel('Time (s)')
    ylabel('Voltage (uV)')
    title(['Baseline --> ', stim_names{i}])
    ylim([ylim_low,ylim_high])
end
sgtitle(['Complex Burst Firing: ',patient_id])

if save_fig == true
    fig_path = strcat(patient_path,'\figures\',f.Name,'_',patient_id);
    savefig(f,fig_path)
end

end 