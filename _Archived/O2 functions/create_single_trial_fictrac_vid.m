function create_single_trial_fictrac_vid(vidSaveDir, imgSaveDir, ftData, sid, tid, varargin)
%=======================================================================================================
% CREATE A MOVIE WITH BEHAVIOR VID COMBINED WITH FICTRAC DATA
%
%
% INPUTS:
%       vidSaveDir  = the directory containing the source vids and optic flow data file
%
%       imgSaveDir  = the directory containing your imaging data and metadata (specifically, requires
%                     files analysisMetadata.mat and ROI_Data_Avg.mat in this directory)
%
%       ftData      = struct containing processed FicTrac data
%
%       sid         = the session ID of the video you want to process.
%
%       tid         = the trial ID of the video you want to process.
%
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       'FrameRate' = (default: 25) the frame rate that the behavior video was recorded at
%
%       'OutputDir' = (default: vidSaveDir) the full path to the directory to save the output video in
%
%       'SmoothWin' = (default: 5) the size of the smoothing window to use on the movement/imaging data
%
%========================================================================================================
try
addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster
addpath(genpath('/home/mjm60/DownloadedFunctions')) % if running on O2 cluster

% Parse optional arguments
p = inputParser;
addParameter(p, 'FrameRate', 25);
addParameter(p, 'OutputDir', vidSaveDir);
addParameter(p, 'SmoothWin', 3);
parse(p, varargin{:});
frameRate = p.Results.FrameRate;
outputDir = p.Results.OutputDir;
smWin = p.Results.SmoothWin;

% Create output diectory if it doesn't already exist
if ~isdir(outputDir)
    mkdir(outputDir)
end

% Load analysis metadata
load(fullfile(imgSaveDir, 'analysisMetadata.mat'), 'analysisMetadata') % --> 'analysisMetadata' struct
vidName = ['sid_', num2str(sid), '_tid_', pad(num2str(tid), 3, 'left', '0')];
vidFile = fullfile(vidSaveDir, [vidName, '.avi']);

% Prepare ficTrac variables
currX = smooth(ftData.mmXYdata(:,1,tid), smWin*3);
currY = smooth(ftData.mmXYdata(:,2,tid), smWin*3);
currFWSpeed = smooth(ftData.mmFWSpeed(:,tid), smWin*3);
% currFWSpeed = smooth(ftData.mmSpeedData(:,tid), smWin*3);
currYawSpeed = smooth(ftData.dHD(:,tid), smWin*3);
currHD = ftData.HD(:,tid);
nROIs = size(ftData.goodFlNorm, 3);

% Calculate shade volumes if applicable
if ~isempty(analysisMetadata.stimOnsetTimes)
    shadeTimes = [analysisMetadata.stimOnsetTimes(tid), ... 
                  analysisMetadata.stimOnsetTimes(tid) + analysisMetadata.stimDurs(tid)];
end

% Calculate frame and volume times
volTimes = (1:analysisMetadata.nVolumes)' ./ analysisMetadata.volumeRate;
frameTimes = (1:analysisMetadata.nFrames)' ./ frameRate;

% Create VidReader and VidWriter
myVid = VideoReader(vidFile);
myVidWriter = VideoWriter(fullfile(vidSaveDir, ['Combined_', vidName, '_With_FicTrac_Plots']), 'Motion JPEG AVI');
myVidWriter.FrameRate = frameRate;
open(myVidWriter)

h = figure(10);

for iFrame = 1:ftData.frameCounts(tid)
    
    currFrame = uint8(readFrame(myVid));
    clf
    
    % Create figure
    screenSize = [1 1 1920 1080];
    h.Color = [1 1 1];
    h.OuterPosition = [50 100 (screenSize(3) - 150) screenSize(4) - 150];
    
    % Create axes
    clf
    M = 0.006;
    axVid = subaxis(2,4, 1, 'S', 0, 'M', M, 'PB', 0.00, 'PR', 0.01);
    axMove = subaxis(2,4, 5, 'S', 0, 'M', M, 'PB', 0.03, 'PR', 0.01);
    axFW = subaxis(2,4, 2:4, 'S', 0, 'M', M, 'PB', 0.06, 'PR', 0.05, 'PL ', 0.03);
    axYaw = subaxis(2,4, 6:8, 'S', 0, 'M', M, 'PB', 0.07, 'PR', 0.05, 'PL ', 0.04);
    
    % Movie frame plot
    imshow(currFrame, [], 'Parent', axVid);
    axVid.XTickLabel = [];
    axVid.Title.String = ['Trial #', num2str(tid), '  -  ', analysisMetadata.trialType{tid}];
    
    % Movement map plot
    axes(axMove);
    hold on
    nSegs = analysisMetadata.trialDuration;
    vectorLen = floor(numel(currX) / nSegs);
    cm = jet(nSegs);
    for iSeg = 1:nSegs
        if iSeg == 1
            padNum = 1;
        else
            padNum = 0;
        end
        currSegInd = (padNum + ((iSeg-1)*vectorLen)):(iSeg * vectorLen);
        plot(axMove, currX(currSegInd), currY(currSegInd), 'Color', cm(iSeg, :), 'LineWidth', 2)
    end
    axis equal;
    lims = 1.1 * max(abs([ylim(axMove), xlim(axMove)]));
    if lims < 1
        lims = 1;
    end
    xlim([-lims lims])
    ylim([-lims lims])
    legend({'2D movement (mm)'}, 'FontSize', 11, 'Autoupdate', 'off', 'location', 'best')
    x = currX(iFrame);
    y = currY(iFrame);
    [arrowVec(1), arrowVec(2)] = pol2cart(deg2rad(currHD(iFrame)), 0.5);
    arrow([x - arrowVec(1)/2, y-arrowVec(2)/2], [x + arrowVec(1)/2, y + arrowVec(2)/2], 'FaceColor', 'green');
    
    % Plot FW move speed + fluorescence for current trial
    axes(axFW); hold on
    axFW.FontSize = 14;
    plotColors = default_plot_colors();
    fwSpeedSmooth = smooth(smooth(currFWSpeed(:), smWin*3), smWin*3);
    yyaxis left
    lgdObj = plot(frameTimes, fwSpeedSmooth, 'LineWidth', 2, 'Color', plotColors(1,:));
    plot(frameTimes(iFrame), fwSpeedSmooth(iFrame), 'd', 'markersize', 10, 'Color', 'k', 'markerfacecolor', 'k');
    ylabel('FW move speed (mm/sec)');
    lgdStr = {'FW move speed'};
    ylim([min(currFWSpeed), max(currFWSpeed)]);
    if ~isempty(analysisMetadata.stimOnsetTimes)
        plot_stim_shading(shadeTimes, 'Axes', axFW);
    end
    yyaxis right
    for iROI = 1:nROIs
        currFl = smooth(ftData.goodFlNorm(:,tid,iROI),smWin);
        lgdObj(end + 1) = plot(volTimes, currFl, ':', 'LineWidth', 3, 'Color', plotColors(iROI + 1, :));
        lgdStr{end + 1} = ['ROI ', num2str(iROI)];
        [~, idx] = min(abs(analysisMetadata.volFrames - iFrame));
        plot(volTimes(idx), currFl(idx), 'p', 'markersize', 14, 'Color', 'k', 'markerfacecolor', 'k');
    end
    ylabel('Raw fluorescence (AU)');
    legend(lgdObj, lgdStr, 'autoupdate', 'off');
    ylim([min(squeeze(ftData.goodFlNorm(:,tid, :))), max(squeeze(ftData.goodFlNorm(:, tid, :)))])
    axFW.FontSize = 14;
    xlim([0 analysisMetadata.trialDuration])
    
    % Plot yaw speed + fluorescence for current trial
    axes(axYaw); hold on
    axYaw.FontSize = 14;
    plotColors = default_plot_colors();
    yawSpeedSmooth = smooth(smooth(currYawSpeed,smWin*3),smWin*3);
    yyaxis left
    lgdObj = plot(frameTimes, yawSpeedSmooth, 'LineWidth', 2, 'Color', plotColors(1,:));
    ylabel('Yaw speed (deg/sec)');
    lgdStr = {'Yaw speed'};
    ylim([min(currYawSpeed(:)), max(currYawSpeed(:))]);
    if ~isempty(analysisMetadata.stimOnsetTimes)
        plot_stim_shading(shadeTimes, 'Axes', axYaw);
    end
    plot(frameTimes(iFrame), yawSpeedSmooth(iFrame), 'd', 'markersize', 10, 'Color', 'k', 'markerfacecolor', 'k');
    yyaxis right
    for iROI = 1:nROIs
        currFl = smooth(ftData.goodFlNorm(:,tid,iROI),smWin);
        lgdObj(end + 1) = plot(volTimes, currFl, ':', 'LineWidth', 3, 'Color', plotColors(iROI + 1, :));
        lgdStr{end + 1} = ['ROI ', num2str(iROI)];
        [~, idx] = min(abs(analysisMetadata.volFrames - iFrame));
        plot(volTimes(idx), currFl(idx), 'p', 'markersize', 14, 'Color', 'k', 'markerfacecolor', 'k');
    end
    ylabel('Raw florescence (AU)');
    legend(lgdObj, lgdStr, 'autoupdate', 'off');
    xlabel('Time (sec)')
    ylim([min(squeeze(ftData.goodFlNorm(:,tid, :))), max(squeeze(ftData.goodFlNorm(:, tid, :)))])
    axYaw.FontSize = 14;
    xlim([0 analysisMetadata.trialDuration])
    
    % Write frame to video(s)
    writeFrame = getframe(h);
    writeVideo(myVidWriter, writeFrame);
end%iFrame

close(myVidWriter)
catch ME
    write_to_log(getReport(ME), mfilename);
end%try
end%function