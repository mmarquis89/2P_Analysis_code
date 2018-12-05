function [goodTrials, frameCounts, totalFrameCount, nFramesTarget] = frame_count_check(parentDir, sid, frameRate, trialDuration)
%===============================================================================================================================
% Checks each individual trial in the experiment for dropped frames in the behavior video, and returns some useful
% information about the frame counts. This function does not count the frames itself (that should already have been done
% during the initial video processing).
%
% INPUTS:
%       parentDir       = the path to the directory containing the .mat files with the frame count information. The two files
%                         should be named as follows: "sid_X_AllTrials_frameCountLog.mat" and "sid_X_frameCountLog.mat".
%
%       sid             = the session ID number of the trials you want to process.
%
%       frameRate       = the frame rate that the behavioral video was acquired at.
%
%       trialDuration   = the total duration of each trial in seconds.
%
% OUTPUTS:
%       goodTrials      = a 1 x nTrials logical vector specifying trials that do not have any dropped frames
%
%       frameCounts     = a 1 x nTrials vector containing the total number of video frames for each trial
%
%       totalFrameCount = the total number of frames across all trials
%
%       nFramesTarget   = the number of frames that each trial should have
%
%===============================================================================================================================

% Load total number of frames in concatenated video
allTrials = load(fullfile(parentDir, ['sid_', num2str(sid), '_AllTrials_frameCountLog.mat']));
totalFrameCount = allTrials.frameCount;

if totalFrameCount > 1
    
    % Load individual trial frame count data
    matFileName = fullfile(parentDir, ['sid_', num2str(sid), '_frameCounts.mat']);
    txtFileName = fullfile(parentDir, ['sid_', num2str(sid), '_frameCounts.txt']);
    if exist(matFileName, 'file')
        individualTrialFrameCounts = load(matFileName);
        frameCounts = [individualTrialFrameCounts.frameCounts.nFrames];
    else
        frameCounts = [];
        myFile = fopen(txtFileName);
        currLine = fgetl(myFile);
        while ischar(currLine)
            tid = str2double(regexp(currLine, '(?<=tid_)...', 'match'));
            frameCounts(tid) = str2double(regexp(currLine, '.*(?=,)', 'match'));
            currLine = fgetl(myFile);
        end
        fclose(myFile);
    end
    
    % Calculate correct frame count
    nFramesTarget = frameRate * trialDuration;
    
    % Check against each trial
    goodTrials = (frameCounts == nFramesTarget);
else
    frameCounts = 0;
    nFramesTarget = 0;
    goodTrials = [];
end


end