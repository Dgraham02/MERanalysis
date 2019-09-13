% Global Function
% Purpose: Extract start/stop times of stimulation from locations of
% stimulation artifacts 

function [start_art, stop_art] = start_stop_times(locs_art)

% Find Start & Stop times of stimulation from artifacts 
diff_stop = diff(locs_art); % find difference between locs in forward direction 
stop1 = diff_stop>2; % find differences greater than 2 
stop2 = [stop1,1]; % add the stop point of last artifact (not captured with diff)
stopB = locs_art.*stop2; % array of zeros and stop points for plotting w/ locs data 
stop_art = locs_art(stopB>0); % locs of stop points ***

diff_start = fliplr(abs(diff(fliplr(locs_art)))); % flip array to find difference in reverse direction for ...
% capturing start points; then flip array back to same direction as locs array 
start1 = diff_start>2; % find differences greater than 2 
start2 = [1,start1]; % add start point of first artifact (not captured w/ diff)
startB = locs_art.*start2; % array of zeros and start point for plotting w/ locs data 
start_art = locs_art(startB>0); % locs of start points ***

end