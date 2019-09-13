% Global Function
% Purpose: Compare Stimulation artifacts to action potentials 
% Notes:

function [f] = fig_art_vs_ap(general_data, post_hp_data, clean_data, save_fig)

% Define variables from structures
patient_path = general_data.patient_path;
patient_id = general_data.patient_id;
Fs = general_data.Fs;
time = general_data.time;
spikeData = post_hp_data.spikeData;
pks_art_remove = clean_data.pks_art_remove;
locs_art_remove = clean_data.locs_art_remove;
pks = clean_data.pks;
locs1 = clean_data.locs1;


f = figure('Name','Art vs AP');clf;
start = 800;
stop = 820;
t1 = (start-time(1))*Fs;
t2 = (stop-time(1))*Fs;
plot(time(t1:t2),spikeData(t1:t2))
hold on
plot(locs_art_remove,pks_art_remove,'o')
plot(locs1,pks,'go')
xlim([start,stop])
xlabel('Time(s)')
ylabel('Voltage(uV)')
title(['Action Potential vs Stimulation Artifact ',patient_id],'Fontsize',16)
legend('Signal', 'Stimulation Artifact', 'Action Potential')

% Removing grey space around image 
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

if save_fig == true
    fig_path = strcat(patient_path,'\figures\',f.Name,'_',patient_id);
    savefig(f,fig_path)
end

end