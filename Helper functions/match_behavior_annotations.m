function behaviorVols = match_behavior_annotations(infoStruct)
%===================================================================================================
%
% Takes trial annotation data and matches it up to imaging volumes (since it is originally annotated 
% by video frames) so it can be used with the GCaMP data. Only necessary argument is infoStruct, 
% but it must contain a number of specific data fields.
%
% INPUTS (infoStruct fields):
%       .nTrials = total number of trials in the experiment.
%
%       .nVolumes = total number of imaging volumes per trial
%
%       .volumeRate = rate of imaging volume acquisition in volumes/sec
%
%       .goodTrials = logical vector indicating trials that aren't missing any behavior video frames
%
%       .trialAnnotations = cell array of length nTrials with a table in each cell containing 
%                           frame-by-frame annotation data for that trial. Table must have fields:
%           .frameTime = trial time of each frame in seconds
%           .actionNums = annotations as numbers of specific actions
%
%       .behaviorLabels = cell array of labels for each unique actionNum value in trialAnnotations
%
% OUTPUTS:
%       behaviorVols = 3D array with dimensions [trial, actionNum, volume]
%
%===================================================================================================

% Pull required data out of infoStruct
nTrials = infoStruct.nTrials;
nVolumes = infoStruct.nVolumes;
volumeRate = infoStruct.volumeRate;
goodTrials = infoStruct.goodTrials;
trialAnnotations = infoStruct.trialAnnotations;
behaviorLabels = infoStruct.behaviorLabels;

behaviorVols = [];
for iTrial = 1:nTrials
    if goodTrials(iTrial)
        
        % Match frame times to volumes
        volTimes = (1:nVolumes)' ./ volumeRate;
        frameTimes = trialAnnotations{find(goodTrials, 1)}.frameTime;
        for iVol = 1:nVolumes
            [~, volFrames(iVol)] = min(abs(frameTimes - volTimes(iVol)));
            volFrames = volFrames';
        end
        
        % Pull out action numbers for each volume
        currBehaviors = trialAnnotations{iTrial}.actionNums;
        volBehaviors = currBehaviors(volFrames);
        
        % Identify volume actions
        for iBehav = 1:length(behaviorLabels)
            behaviorVols(iTrial, iBehav, :) = (volBehaviors == (iBehav - 1)); % --> [trial, behavior, volume]
        end
    else
        % So data from invalid trials won't ever be matched to a behavior state
        behaviorVols(iTrial, :, :) = 0;
    end
end

end