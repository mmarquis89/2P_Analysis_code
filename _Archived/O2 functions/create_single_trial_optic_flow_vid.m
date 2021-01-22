function create_single_trial_optic_flow_vid(vidSaveDir, sid, tid, frameCountFile, varargin)
%=======================================================================================================
% CREATE A MOVIE WITH BEHAVIOR VID COMBINED WITH OPTIC FLOW DATA
%
% Uses previously defined behavior vid ROI to create a movie that combines the behavior video with a
% sliding plot of optic flow around the fly to make behavioral annotation in ANVIL easier. Requires
% a .txt file with comma-separated frame count data for all the trials in <frame count>,<sid_X_tid_XXX>
% format.
%
% INPUTS:
%       vidSaveDir   = the directory containing the source vids and optic flow data file
%
%       sid         = the session ID of the video you want to process.
%
%       tid         = the trial ID of the video you want to process.
%
%       frameCountFile =  .txt file with frame counts for each trial in [frames],[tid] format
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       'FrameRate' = (default: 25) the frame rate that the behavior video was recorded at
%
%       'OutputDir' = (default: vidSaveDir) the full path to the directory to save the output video in
%
%       'FlowDataDir' = (default: vidSaveDir) the full path to the directory containing optic flow data
%
%
%========================================================================================================
try
addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

% Parse optional arguments
p = inputParser;
addParameter(p, 'FrameRate', 25);
addParameter(p, 'OutputDir', vidSaveDir);
addParameter(p, 'FlowDataDir', vidSaveDir);
parse(p, varargin{:});
frameRate = p.Results.FrameRate;
outputDir = p.Results.OutputDir;
flowDataDir = p.Results.FlowDataDir;

disp(outputDir)

% Create output directory if necessary
if ~isdir(outputDir)
    mkdir(outputDir)
end

% Load optic flow data
load(fullfile(flowDataDir, ['sid_', num2str(sid), '_flow_data_norm.mat'])); % 'flyFlowNorm'

% Load frame count log 
frameCounts = load_frame_counts(vidSaveDir, frameCountFile, sid);

% Figure out how many trials will need to be plotted at some point in this video
preAlignTime = 4;
postAlignTime = 16;
xLimFrames(1) = round(preAlignTime * frameRate);
xLimFrames(2) = round(postAlignTime * frameRate);
[preTrialCount, postTrialCount, preAlignFrames, postAlignFrames] = ...
get_required_trials(preAlignTime, postAlignTime, frameCounts, frameRate, tid);

% Get optic flow data from current and any other required trials
currTrialFlow = flyFlowNorm{tid};
nFrames = numel(currTrialFlow);
preTrialFlow = [];
if preTrialCount > 0 
    
    for iTrial = 1:preTrialCount
        currTrial = tid - iTrial;
        preTrialFlow = [flyFlowNorm{currTrial}, preTrialFlow];
    end
end
postTrialFlow = [];
if postTrialCount > 0 
    for iTrial = 1:postTrialCount
        currTrial = tid + iTrial;
        postTrialFlow = [postTrialFlow, flyFlowNorm{currTrial}];
    end
end

% Pull out the relevant flow data frames for other trials
if ~isempty(preTrialFlow)
    preTrialFlow = preTrialFlow(end-preAlignFrames+1:end);
end
if ~isempty(postTrialFlow)
    postTrialFlow = postTrialFlow(1:postAlignFrames); 
end
plotFlowData = [preTrialFlow, currTrialFlow, postTrialFlow];

% Calculate locations of all trial boundaries that need to be plotted
trialBounds  = [preAlignFrames + 1, preAlignFrames + nFrames + 1]; % Current trial bounds
if preTrialCount > 1
    for iTrial = 1:preTrialCount - 1
        trialBounds(end + 1) = (preAlignFrames + 1) - sum(frameCounts(tid-iTrial:tid-1));
    end
end
if postTrialCount > 1
    for iTrial = 1:postTrialCount - 1
        trialBounds(end + 1) = (preAlignFrames + nFrames + 1) + sum(frameCounts(tid+1:tid+iTrial));
    end
end

% Load current trial's video
currVid = ['sid_', num2str(sid), '_tid_', pad(num2str(tid), 3, 'left', '0'), '.avi'];
trialVid = VideoReader(fullfile(vidSaveDir, currVid));
trialFrames = [];
frameCount = 1;
while hasFrame(trialVid)
    currFrame = readFrame(trialVid);
    trialFrames(:, :, frameCount) = uint8(currFrame(:,:,1));
    frameCount = frameCount + 1;
end

% Create vidWriter
outputFileName = ['sid_', num2str(sid), '_tid_', pad(num2str(tid), 3, 'left', '0'), '_With_Optic_Flow'];
myVidWriter = VideoWriter(fullfile(outputDir, outputFileName), 'Motion JPEG AVI');
myVidWriter.FrameRate = frameRate;
open(myVidWriter)

% Plot and write each frame
for iFrame = (preAlignFrames + 1):(preAlignFrames + nFrames)
    
    disp(num2str(iFrame))
    currVidFrame = trialFrames(:,:, iFrame-preAlignFrames);
    xSize = size(currVidFrame, 2);
    
    % Create figure
    screenSize = [1 1 1824 1026];
    h = figure(10);clf
    h.OuterPosition = [50 50 (xSize) screenSize(4)-50];
    
    % Movie frame plot
    ax = axes('Units', 'Normalized', 'Position', [0 0 1 0.7]);
    imshow(currVidFrame, []);
    axis normal;
    ax.Units = 'Pixels';
    minPos = min(ax.Position(3:4));
    ax.Position(3:4) = [minPos, minPos];
    h.Position(3) = ax.Position(3);
    axis off
    set(gca, 'xticklabel', []);
    
    % Calculate optic flow xLims
    xL = [iFrame - xLimFrames(1), iFrame + xLimFrames(2)];
    if xL(1) <= 0
        xL(1) = 0;
        xL(2) = preAlignFrames + postAlignFrames;
    elseif xL(2) > numel(plotFlowData)
        xL(1) = numel(plotFlowData) - (preAlignFrames + postAlignFrames);
        xL(2) = numel(plotFlowData);
    end
    
    % Optic flow plot
    axes('Units', 'Pixels', 'Position', [0 ax.Position(4) ax.Position(3) (h.Position(4) - ax.Position(4))]);
    hold on
    plot(smoothdata(plotFlowData(2:end), 'gaussian', 5), 'color', 'm');    % Plot fly movmement ROI flow
    ylim([0 1]);
    plot([iFrame, iFrame], ylim(), 'LineWidth', 2, 'color', 'r')
    for iTrial = 1:length(trialBounds)
        plot([trialBounds(iTrial), trialBounds(iTrial)], ylim(), 'LineWidth', 2, 'color', 'k')
    end
    set(gca, 'xticklabel', []);
    xlim(xL);
    ylabel('Optic flow (au)');
    lgd = legend('Fly movement');
    lgd.LineWidth = 1;
    
    drawnow()
    
    % Write frame to video
    writeFrame = getframe(h);
    if iFrame > (preAlignFrames + 1)
        % Something's messed up with the first frame's size...
        writeVideo(myVidWriter, writeFrame);
    end
    if iFrame == (preAlignFrames + 2)
        % ...so I'm replacing it with another copy of the second frame
        writeVideo(myVidWriter, writeFrame);
    end
    
end%iFrame

close(myVidWriter)

catch ME
    write_to_log(getReport(ME), mfilename);
end%try
end%function



% ===============================================================================================
% LOCAL FUNCTIONS
% ===============================================================================================


%-------------------------------------------------------------------------------------
% Load frame count log and extract frame counts for each trial in the current session
%-------------------------------------------------------------------------------------
function frameCounts = load_frame_counts(vidSaveDir, frameCountFile, sid)

frameCountFile = fopen(fullfile(vidSaveDir, frameCountFile), 'r');
frameCounts = [];
currLine = fgetl(frameCountFile);
while ischar(currLine)
    currSid = str2double(regexp(currLine, '(?<=sid_).', 'match'));
    currTid = str2double(regexp(currLine, '(?<=tid_)...', 'match'));
    if currSid == sid
        frameCounts(currTid) = str2double(regexp(currLine, '.*(?=,)', 'match'));
    end
    currLine = fgetl(frameCountFile);
end
fclose(frameCountFile);

end 

%------------------------------------------------------------------------------------
% Figure out how many trials will need to be plotted at some point in this video
% based on the number of frames that will need to be on screen at any one time and 
% the number of frames of each individual trial video.
%------------------------------------------------------------------------------------
function [preTrialCount, postTrialCount, preAlignFrames, postAlignFrames] = ...
         get_required_trials(preAlignTime, postAlignTime, frameCounts, frameRate, tid)

preAlignFrames = round(preAlignTime * frameRate);
postAlignFrames = round(postAlignTime * frameRate);
preRunningCount = 0;
preTrialCount = 0;
while preRunningCount < preAlignFrames
    if (tid - preTrialCount) > 1
        preRunningCount = preRunningCount + frameCounts(tid - preTrialCount - 1);
        preTrialCount = preTrialCount + 1;
    else
        preAlignFrames = preRunningCount;
    end
end
postRunningCount = 0;
postTrialCount = 0;
while postRunningCount < postAlignFrames
    if (tid + postTrialCount) < numel(frameCounts)
        postRunningCount = postRunningCount + frameCounts(tid + postTrialCount + 1);
        postTrialCount = postTrialCount + 1;
    else
        postAlignFrames = postRunningCount;
    end
end
end
