% Global Function
% Purpose: Separate Data Based on Stim On/Off
% Notes: 
function [OnOff] = separate_on_off(general_data,post_hp_data,cleanSpikeData)
time = general_data.time;
Fs = general_data.Fs;
start_art = post_hp_data.start_art;
stop_art = post_hp_data.stop_art;
blocks = post_hp_data.blocks;
numBlocks = post_hp_data.numBlocks;

start_art2 = round((start_art-time(1))*Fs); % convert start_art from time to samples and normalize to current block 
stop_art2 = round((stop_art-time(1))*Fs); % convert stop_art from time to samples and normalize to current block

Data = cleanSpikeData; % data to be used --> filtered/stim artifacts removed 
%Data = sig; % use unfiltered data 

% OnOff('block', 'stim#', 'On/Off','sig/time') On=1,Off=2
OnOff = cell(numBlocks,(length(start_art)+2),2,2); % initialize cell array to number of stims +2 baselines 
for ii = 1:length(blocks) 
    for i = 1:length(start_art(start_art<blocks{ii,2}(end)))-(ii-1)*10 %adjust stim trial length per block
        if start_art(i)<=blocks{ii,2}(end)
            scaling_index = i + ((ii-1)*10); % scale start_art index to each block 
            OnOff{ii,i,1,1} = Data(start_art2(scaling_index):stop_art2(scaling_index)); % stim ON spikeData (sig)
            OnOff{ii,i,1,2} = time(start_art2(scaling_index):stop_art2(scaling_index)); % stim ON spikeData (time)
            if i < (length(start_art(start_art<blocks{ii,2}(end)))-(ii-1)*10) % do not attempt to index stim 11 when there are only 10 
                OnOff{ii,i,2,1} = Data(stop_art2(scaling_index):start_art2(scaling_index+1));% stim OFF spikeData (sig)
                OnOff{ii,i,2,2} = time(stop_art2(scaling_index):start_art2(scaling_index+1));% stim OFF spikeData(time)
            else
            
                block_start_index = round((blocks{ii,2}(1)-time(1))*Fs); % first index of block 
                block_stop_index = round((blocks{ii,2}(end)-time(1))*Fs); % last index of block
                OnOff{ii,i+1,2,1} = Data(block_start_index+1:start_art2(scaling_index-(i-1))); % pre-stim baseline (sig)
                OnOff{ii,i+2,2,1} = Data(stop_art2(scaling_index): block_stop_index); % post-stim baseline (sig)
                OnOff{ii,i+1,2,2} = time(block_start_index+1:start_art2(scaling_index-(i-1))); % pre-stim baseline (time)
                OnOff{ii,i+2,2,2} = time(stop_art2(scaling_index):block_stop_index); % post-stim baseline (time)
                
            end
        else
        end
    end
end
