% Global Function
% Purpose: Determine a quantiative threshold for Action Potentials
% Notes: 
function [APtime, pks, locs1] = find_APs(general_data, post_hp_data, threshold, Data)

time = general_data.time;
start_art = post_hp_data.start_art;

figure(5); clf;
hold on
plot(time,Data)
thresh_line = refline(0,threshold); % plot threshold line 
thresh_line.Color = 'r';
for i = 1:length(start_art)
    line([start_art(i) start_art(i)], [min(Data)*1.2,max(Data)*1.2],'Color','green')
end
xlabel('Time(s)')
ylabel('Voltage(uV)')
title('Threshold & Ramp')
legend('Filtered Data','Threshold Line','Start Times of Stimulation','Location','southeast')

% Find Peaks 
APtime = 2/1000; % duration of an action potential = 2ms
[pks, locs1] = findpeaks(Data,time,'MinPeakHeight',threshold,'MinPeakDistance',APtime);
