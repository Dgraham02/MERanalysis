% Global Function: Separate data into recording blocks using baselineTimes

function [blocks, numBlocks] = separate_blocks(general_data, startStop)

C1 = general_data.C1;
time = general_data.time;
baselineTimes = general_data.baselineTimes;
Fs = general_data.Fs;
tRaw = general_data.tRaw;

numBlocks = length(startStop(:,1)); % number of rows in 'startStop'

uVC1 = C1*1e6; % convert signal from V to uV 

blocks = cell(numBlocks,2); % initialize cell array
for i = 1:numBlocks % number of blocks = 4 
    startIndex = round(baselineTimes(startStop(i,1))*Fs); % index of 'start' baseline marker
    stopIndex = round(baselineTimes(startStop(i,2))*Fs); % index of 'stop' baseline marker 
    
    if i==numBlocks % capture last baseline recording after stim 
        stopIndex = round(time(end)*Fs);
    end
    blocks{i,1} = uVC1(startIndex:stopIndex); % assign voltage values to column 1 of 'blocks' 
    blocks{i,2} = tRaw(startIndex:stopIndex);% assign time values to column 2 of 'blocks'
end
end