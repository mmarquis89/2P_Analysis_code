function make_optic_flow_vid(parentDir, inputVid, FRAME_RATE, varargin)
%=======================================================================================================
% CREATE A MOVIE WITH BEHAVIOR VID COMBINED WITH OPTIC FLOW DATA
%
% Uses previously defined behavior vid ROI to create a movie that combines the behavior video with a
% sliding plot of optic flow around the fly to make behavioral annotation in ANVIL easier. 
%
% INPUTS:
%       parentDir   = the directory containing the optic flow video
%
%       inputVid    = the name of the input video file to use (minus the .mp4 extension)
%
%       FRAME_RATE  = the frame rate of the behavior video
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       'ROIfile'   = (default: []) the full path to the .mat file containing the ROI data mask. Can use [] 
%                     to have the function prompt the user to select a file.
% 
%       'OutputDir' = (default: parentDir) the full path to the directory to save the output video in
%
%       'FrameCountFile' = (default: []) the full path to the .mat file containing frame count data
%
%       'OpticFlowFile' = (default: []) the full path to a .mat file containing optic flow data
%
%       'OutputFileName' = (default: [inputVid, '_With_Optic_Flow']) name for the output video file
%
%========================================================================================================

addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

% Parse optional arguments
p = inputParser;
addParameter(p, 'ROIfile', [])
addParameter(p, 'OutputDir', parentDir);
addParameter(p, 'FrameCountFile', []);
addParameter(p, 'OpticFlowFile', []);
addParameter(p, 'OutputFileName', [inputVid, '_With_Optic_Flow']);
parse(p, varargin{:});
outputDir = p.Results.OutputDir;
roiDataFile = p.Results.ROIfile;
frameCountFile = p.Results.frameCountFile;
flowDataFile = p.Results.OpticFlowFile;
outputFileName = p.Results.OutputFileName;

% Load ROI data file
if isempty(roiDataFile)
    roiDataFile = uigetfile([parentDir, '\*.mat'], 'Select an ROI data file');
end
load(roiDataFile);

% Load frame count log
if isempty(frameCountFile)
    frameCountFile = uigetfile([parentDir, '\*.mat'], 'Select a frame count file');
end
individualVidFrameCounts = load(frameCountFile);
frameCounts = [];
for iTrial = 1:length(individualVidFrameCounts.frameCounts)
    if isempty(individualVidFrameCounts.frameCounts(iTrial).nFrames)
        frameCounts(iTrial) = 0;
    else
        frameCounts(iTrial) = individualVidFrameCounts.frameCounts(iTrial).nFrames;
    end   
end
framesPerTrial = mode(frameCounts);

% Calculate optic flow for each movie frame (unless file already exists)    
vidFile = fullfile(parentDir, [inputVid, '.mp4']);
 if isempty(dir(flowDataFile))
     
    myVid = VideoReader(vidFile);
    frameCount = 0;
    opticFlow = opticalFlowFarneback;
    while hasFrame(myVid)
        
        frameCount = frameCount + 1;        
        currFrame = readFrame(myVid);
        currFrame = currFrame(:,:,1);
        
        % Calculate optic flow within each ROI
        currFrameFlowData = estimateFlow(opticFlow, currFrame);
        currFlowFly = currFrameFlowData.Magnitude;
        currFlowFly(~roiData) = nan;
        meanFlowMag(frameCount) = nanmean(nanmean(currFlowFly));
        
    end
    nFrames = frameCount;    
    
    % Calculate frame times
    frameTimes = (1:1:nFrames) ./ FRAME_RATE;
    
    % Normalize flow magnitudes
    flyFlow = meanFlowMag;
    flyFlowNorm = flyFlow ./ max(flyFlow(2:end)); % First frame is artificially high so don't use that
    
    % Save optic flow data
    savefast(fullfile(outputDir, ['sid_', num2str(sid), '_optic_flow_data.mat']), 'flyFlowNorm', 'nFrames', 'frameTimes');
 
 else
     % Just load the existing data if it exists
    load(flowDataFile);
 end
 
% Recreate vidReader
myVid = VideoReader(vidFile);

% Create vidWriter
myVidWriter = VideoWriter(fullfile(outputDir, outputFileName), 'MPEG-4');
myVidWriter.FrameRate = FRAME_RATE;
open(myVidWriter)

% Plot and write each frame
for iFrame = 1:nFrames
    
    if ~mod(iFrame, framesPerTrial)
        disp(['Plotting trial #' num2str(iFrame/framesPerTrial)])
    end
    currFrame = uint8(readFrame(myVid));
    ySize = size(currFrame, 1);
    xSize = size(currFrame, 2);
    
    % Create figure
    screenSize = [1 1 1824 1026];
    h = figure(10);clf
    h.OuterPosition = [50 50 (xSize) screenSize(4)-50];
    
    % Movie frame plot
    ax = axes('Units', 'Normalized', 'Position', [0 0 1 0.7]);
    imshow(currFrame, []);
    axis normal;
    ax.Units = 'Pixels';
    minPos = min(ax.Position(3:4));
    ax.Position(3:4) = [minPos, minPos];
    h.Position(3) = ax.Position(3);
    axis off
    set(gca, 'xticklabel', []);    
    
    % Calculate optic flow xLims
    currFrameTime = iFrame * (1/FRAME_RATE);
    preFrameTime = 4;
    postFrameTime = 16;
    xL = [currFrameTime - preFrameTime, currFrameTime + postFrameTime];
    if xL(1) <= 0
        xL(1) = 0;
        xL(2) = preFrameTime + postFrameTime;
    elseif xL(2) > frameTimes(end)
        xL(1) = frameTimes(end) - (preFrameTime + postFrameTime);
        xL(2) =  frameTimes(end);
    end
    
    % Optic flow plot
    trialBoundTimes = [];
    runningCount = 1;
    for iTrial = 1:(length(frameCounts)-1)
        if runningCount > length(frameTimes) % This prevents an error if the first one or more trials has zero frames
           runningCount = length(frameTimes); 
        end
        trialBoundTimes(iTrial) = frameTimes(runningCount);
        runningCount = runningCount + frameCounts(iTrial);
    end
    trialBoundTimes(end+1) = frameTimes(end);
    
    axes('Units', 'Pixels', 'Position', [0 ax.Position(4) ax.Position(3) (h.Position(4) - ax.Position(4))]);
    hold on
    plot(frameTimes(2:end), smooth(flyFlowNorm(2:end), 5), 'color', 'm');    % Plot fly movmement ROI flow
    plot([currFrameTime, currFrameTime], ylim(), 'LineWidth', 2, 'color', 'r')
    for iTrial = 1:length(trialBoundTimes)
       plot([trialBoundTimes(iTrial), trialBoundTimes(iTrial)], ylim(), 'LineWidth', 2, 'color', 'k') 
    end
    set(gca, 'xticklabel', []);
    xlim(xL);
    ylabel('Optic flow (au)');
    lgd = legend('Fly movmement');
    lgd.LineWidth = 1;
    
    drawnow()
    % Write frame to video
    writeFrame = getframe(h);
    writeVideo(myVidWriter, writeFrame);   
    
end

close(myVidWriter)

end