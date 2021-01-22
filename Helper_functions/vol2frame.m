function frameNums  = vol2frame(volNums, volTimes, frameTimes)
%===================================================================================================
% 
% Converts a numeric array of one or more imaging volumes into video frame numbers.
% 
% INPUTS: 
%       volNums     = the input vector of the volumes that you want to convert to frame numbers 
%       volTimes    = a vector of length nVolumes containing the trial time (in sec) for each volume
%       frameTimes  = same as volTimes, only for behavior video frames
%
% OUTPUTS:
%       frameNums  = the numbers of the frame that is associated with each volume. Note that because 
%                    FRAME_RATE > volumeRate, there are multiple frames that are closest to any 
%                    given volume number, and this function will just return the frame number with 
%                    the frameTimes value closest to the volTimes value. This means that 
%                    vol2frame(frame2vol(X)) ~= X in all cases.
%
%===================================================================================================
frameNums = [];
for iVol = 1:numel(volNums)
    [~, frameNums(iVol)] = min(abs(frameTimes - volTimes(volNums(iVol))));
end

end