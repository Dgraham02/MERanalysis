% Global Function
% Purpose: Validate Start and Stop times derived from stim artifacts 
% Notes: 

function [f] = fig_validate_start_stop(general_data, post_hp_data, save_fig)

% Define variables from structures 
patient_path = general_data.patient_path;
patient_id = general_data.patient_id;
time = general_data.time;
spikeData = post_hp_data.spikeData;
threshold_art = post_hp_data.threshold_art;
locs_art = post_hp_data.locs_art;
pks_art = post_hp_data.pks_art;
start_art = post_hp_data.start_art;
stop_art = post_hp_data.stop_art;


f = figure('Name','Validate Start Stop'); clf;
ax1 = subplot(2,1,1);
plot(time,spikeData,locs_art,pks_art,'o')
thresh_line = refline(0,threshold_art); % plot threshold line 
thresh_line.Color = 'r';

xlim([time(1),time(end)])
xlabel('Time(s)')
ylabel('Voltage(uV)')
sgtitle('Validate Start/Stop Times of Stimulation')

ax2 = subplot(2,1,2);
scale = 500;
for i = 1:length(start_art)
    line([start_art(i) start_art(i)], [-scale,scale],'Color','green')
    line([stop_art(i) stop_art(i)], [-scale,scale],'Color','red')
    patch([start_art(i), stop_art(i), stop_art(i), start_art(i)], [-scale, -scale, scale, scale], 1, 'facecolor', 'r', 'edgecolor', 'none', 'facealpha', .15  );
end
linkaxes([ax1,ax2],'x')

if save_fig == true
    fig_path = strcat(patient_path,'\figures\',f.Name,'_',patient_id);
    savefig(f,fig_path)
end

end