function [trimmedData, preStimData, stimData] = extract_target_volumes(dataArr, onsetTime, volumeRate, baselineDurSec, respDurSec)
%==============================================================================================================
% 
% Extract imaging data from a specified period, and return that and/or just the
% pre-stimulus baseline period or the period during the stimulus.
%
% INPUTS:
%       dataArr        = 5-D imaging data array with dimensions [y, x, plane, volume, trial]
%
%       onsetTime      = trial time in seconds that you want to align the analysis to 
%
%       volumeRate     = the rate of imaging data acquisition in volumes/sec
%                        
%       baselineDurSec = duration of the baseline period in seconds
%
%       respDurSec     = duration of the stimulus/response period in seconds
%
% OUTPUTS:
%       trimmedData    = session data cropped to include just the specified time window
%
%       preStimData    = just the pre-stimulus baseline
% 
%       stimData       = just the post-baseline period
% 
%==============================================================================================================

% Calculate stim period volumes
onsetVol = ceil(onsetTime * volumeRate);
baselineStartVol = onsetVol - floor(baselineDurSec * volumeRate);
respEndVol = floor(onsetVol + (respDurSec * volumeRate));

% Extract data
trimmedData = dataArr(:,:,:, baselineStartVol:respEndVol);
preStimData = dataArr(:,:,:, baselineStartVol:onsetVol-1);
stimData = dataArr(:,:,:, onsetVol:respEndVol);

end