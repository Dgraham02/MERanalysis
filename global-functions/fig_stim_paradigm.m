% Global Function
% Purpose: Show plot of stimulation/recording paradigm (for poster)
% Notes: Start/stop designed for patient-2018-02

function [f] = fig_stim_paradigm(general_data, post_hp_data, save_fig)

% Define variables from structures
patient_path = general_data.patient_path;
patient_id = general_data.patient_id;
Fs = general_data.Fs;
time = general_data.time;
spikeData = post_hp_data.spikeData;
start_art = post_hp_data.start_art;
stop_art = post_hp_data.stop_art; 

ylim_low = min(spikeData); % upper voltage bound 
ylim_high = max(spikeData); % lower voltage bound 

f = figure('Name','Stim Paradigm');clf;
start = 600;
stop = 1080;
t1 = (start-time(1))*Fs;
t2 = (stop-time(1))*Fs;
scale = ylim_high;

plot(time(t1:t2),spikeData(t1:t2))
xlim([start,stop])
ylim([ylim_low,ylim_high+100])
xlabel('Time(s)')
ylabel('Voltage(uV)')
title(['Stimulation/Recording  ', patient_id], 'Fontsize',16)

for i = 1:length(start_art)
    patch([start_art(i), stop_art(i), stop_art(i), start_art(i)], [-scale/2.5, -scale/2.5, scale, scale], 1, 'facecolor', 'r', 'edgecolor', 'none', 'facealpha', .15  );
hold on
end

ydim = 0.6;
dim = [.44 ydim .1 .3];
str = 'Stimulation';
annotation('textbox',dim,'String',str,'FitBoxToText','on');

dim = [.16 ydim .1 .3];
str = 'Baseline';
annotation('textbox',dim,'String',str,'FitBoxToText','on');

dim = [.75 ydim .1 .3];
str = 'Baseline';
annotation('textbox',dim,'String',str,'FitBoxToText','on');

% Removing grey space around image 
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];

if save_fig == true
    fig_path = strcat(patient_path,'\figures\',f.Name,'_',patient_id);
    savefig(f,fig_path)
end

end