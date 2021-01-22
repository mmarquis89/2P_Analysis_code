function trialAnnotations = process_anvil_annotations(sid, parentDir, saveDir, annotationFileName, frameRate, trialDuration)
%========================================================================================================================= 
% READS AND PROCESSES A FRAME-BY-FRAME BEHAVIORAL ANNOTATION FILE CREATED IN ANVIL
%
% This function requires information about the frame counts for both the individual and concatenated video files for 
% the selected session, saved in the parent directory as 'sid_0_frameCountLog.mat' and 
% 'sid_0_AllTrials_frameCountLog.mat', respectively (with the appropriate session ID #).
%
% Processed annotation data will be saved as a .mat file in the specified directory containing a table with the 
% annotation data as well as a cell array with labels for the behavioral codes and a logical vector of valid trials.
%
% INPUTS:
%       sid  =  session ID of the data you want to process.
%       parentDir  =  path to the directory containing the annotation file and the frame count logs.
%       saveDir  =  directory to save processed annotation data in.
%       annotationFileName = name of the Anvil-exported .txt file with annotations.
%       frameRate = frame rate that the behavioral video was acquired at.
%       trial Duration = duration of each trial in seconds (must be the same for all trials).
%
% OUTPUTS:
%       trialAnnotations = a 1 x n cell array (where n is the number of trials in the session), with each cell containing an 
%                          m x 3 table (where m is the number of video frames for that trial) with the following column 
%                          names: [frameNum, actionNums, frameTime]. [actions] contains a behavioral code for each frame, 
%                          corresponding to an entry in the "behaviorLabels" array. [frameTime] is the trial time in 
%                          seconds corresponding to each frame. If there are two tracks in the annotation data, this table
%                          will also have a column with ball stopping annotations: 
%                          [frameNum, actionNums, ballStopNums, frameTime]
%
%       behaviorLabels   = an array of strings containing the various strings that correspond to numbers in the [actionNums]
%                          field of the trialAnnotations tables (0 = first entry in behaviorLabels, etc.) Note hardcoded
%                          values for this variable below. Also saves ballStopLabels if there are two annotation tracks.
%
%       goodTrials       = a 1 x n logical array (where n is the number of trials in the session) indicating which trials 
%                          are missing one or more video frames.
%                           
%==========================================================================================================================

% Remember to update this if annotation coding changes
behaviorLabels = {'None', 'AbdominalContraction', 'Locomotion', 'Grooming', 'IsolatedLegMovement'};
ballStopLabels = {'None', 'WasherMoving', 'BallStopped'};

% Read and parse annotation table data
annotationData = readtable(fullfile(parentDir, annotationFileName));
frameNum = annotationData.Frame;
frameNum = frameNum + 1; % Use 1-indexing for frame numbers
 
% % Remove any extra columns
% annotationData(:,isnan(annotationData{1,:})) = [];

% Have to specifically remove the column called "var7" due to a 2017b change to readtable()
% annotationData.Var7 = [];

% Determine how many tracks there are
nTracks = (size(annotationData, 2) - 2) / 2;

% Process annotation data for each track
actionNums = annotationData.ActionTypes_ActionTypes;
actionNums(actionNums == -1000) = 0;
annotationTable = table(frameNum, actionNums);
ballStopTable = [];
if nTracks == 2   ballStopNums = annotationData.BallStopping_BallStopping;
   ballStopNums(ballStopNums == -1000) = 0;
   ballStopTable = table(ballStopNums);
   annotationTable = [annotationTable, ballStopTable];
end

% Check frame counts for each trial
[goodTrials, frameCounts, allTrialsFrameCount, ~] = frame_count_check(parentDir, sid, frameRate, trialDuration);

% Calculate frame times in seconds for good trials
frameTimes = (1:frameRate*trialDuration) * (1/frameRate);

% Make sure frame counts are consistent with each other and annotation data
assert(sum(frameCounts) == allTrialsFrameCount, 'Error: sum of individual frame counts is not equal to concatenated video frame count');
assert(allTrialsFrameCount == length(frameNum), 'Error: frame number mismatch between annotation data and video file data');

% Separate annotation data into individual trials
nTrials = numel(goodTrials);
trialAnnotations = cell(1,nTrials);
currFrame = 1;
for iTrial = 1:nTrials
    
    lastFrame = currFrame+frameCounts(iTrial)-1;

    % Create column of times to add to table(s)
    if goodTrials(iTrial)
        frameTime = frameTimes';
    else
        frameTime = zeros(1, frameCounts(iTrial))'; % Replace times with zeros for trials without enough frames
    end
    trialAnnotations{iTrial} = annotationTable(currFrame:lastFrame,:);
    trialAnnotations{iTrial} = [trialAnnotations{iTrial}, table(frameTime)];
    currFrame = currFrame + frameCounts(iTrial);
end

% Make sure output data file doesn't already exist for this session, then save processed annotation
% data along with behavior labels
saveFilePath = fullfile(saveDir, ['Annotations.mat']);
% assert(exist(saveFilePath, 'file') == 0, 'Error: a file with this name already exists in the save directory')
if nTracks == 2
    save(saveFilePath, 'trialAnnotations', 'behaviorLabels', 'ballStopLabels', 'goodTrials', '-v7.3');
else
    save(saveFilePath, 'trialAnnotations', 'behaviorLabels', 'goodTrials', '-v7.3');
end
end%function

