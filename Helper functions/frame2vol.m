function volumeNums  = frame2vol(frameNums, volTimes, frameTimes)
%===================================================================================================
% 
% Converts a numeric array of one or more video frame numbers into imaging volumes.
% 
% INPUTS: 
%       frameNums   = the input vector of frame numbers that you want to convert to volumes
%       volTimes    = a vector of length nVolumes containing the trial time (in sec) for each volume
%       frameTimes  = same as volTimes, only for behavior video frames
%
% OUTPUTS:
%       volumeNums  = the numbers of the imaging volumes that most closely match each frameNum
%
%===================================================================================================
volumeNums = [];
for iFrame = 1:numel(frameNums)
    [~, volumeNums(iFrame)] = min(abs(volTimes - frameTimes(frameNums(iFrame))));
end

end