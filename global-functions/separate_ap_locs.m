% Global Function
% Purpose: Create a structure of AP locations separated by On/Off and Block
% Notes: Not catching last stim trial, need to explorer solutions 

function [LFPstruc, locs2] = separate_ap_locs(general_data, ap_data, post_hp_data, OnOff_lfp)

time = general_data.time;
blocks = post_hp_data.blocks;
locs = ap_data.locs;
start_art = post_hp_data.start_art;
stop_art = post_hp_data.stop_art;

% OnOff('block', 'stim#', 'On/Off','sig/time') On=1,Off=2

locs_lfp=NaN; % initialize AP locs to be used for lfp alignment 
locs2 = locs+time(1); % adjust locs to align with correct time 
for on_off = 1:2
    for ii = 1:length(blocks)
        %for i = 1:length(start_art(start_art<blocks{ii,2}(end)))-(ii-1)*10
        for i = 1:9
            place_holder = locs2(locs2 > OnOff_lfp{ii,i,on_off,2}(1) & locs2<OnOff_lfp{ii,i,on_off,2}(end));
            locs_lfp = [locs_lfp, place_holder];
        end
        
        if on_off ==1
            myfield = strcat('on_block',num2str(ii));
        else
            myfield = strcat('off_block',num2str(ii));
        end
        LFPstruc.(myfield) = locs_lfp;
        locs_lfp=NaN;
    end
end

% Save baselines in structure 
for on_off = 2
    for ii = 1:length(blocks)
        for i = 11
            place_holder = locs2(locs2 > OnOff_lfp{ii,i,on_off,2}(1) & locs2<OnOff_lfp{ii,i,on_off,2}(end));
            locs_lfp = [locs_lfp, place_holder];
        end
        
        myfield = strcat('baseline_block',num2str(ii));

        LFPstruc.(myfield) = locs_lfp;
        locs_lfp=1000;
    end
end

figure(10);clf;
block_names = fieldnames(LFPstruc);
for ii = 1:4 % Plot locs in ON stim periods 
    plot(LFPstruc.(block_names{ii}),ones(1,length(LFPstruc.(block_names{ii}))),'ro')
    hold on
end

for ii = 5:8 % Plot locs in OFF stim periods
    plot(LFPstruc.(block_names{ii}),ones(1,length(LFPstruc.(block_names{ii}))),'b*')
end

for ii = 9:12 % Plot baseline periods
    plot(LFPstruc.(block_names{ii}),ones(1,length(LFPstruc.(block_names{ii}))),'b*')
end

scale = 2;
for i = 1:length(start_art) % plot stim patches
    patch([start_art(i), stop_art(i), stop_art(i), start_art(i)], [0, 0, scale, scale], 1, 'facecolor', 'r', 'edgecolor', 'none', 'facealpha', .15  );
end
xlabel('Time(s)')
ylabel('Firing Rate')
end