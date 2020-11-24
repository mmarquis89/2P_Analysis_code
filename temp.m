
%% Group data from several compatible trials into a single block


expID = '20201120-2';
trialNums = [1:7];

flowSmWin = 7;
moveThresh = 0.08;

flType = 'expDff';

sectionBreakTrials = [3 6];
sectionColorLabels = {'5 mM ATP (60 sec pulse)', '5 mM ATP (60 sec pulse)'};%{'100 uM DA', '250 uM DA', '500 uM DA', 'No visual stim', 'Back to visual stim'};
sectionBreakColors = [1:numel(sectionBreakTrials)];
fullColorList = [rgb('magenta'); rgb('green'); rgb('gold'); rgb('cyan') * 0.8; rgb('red') * 0.8];
sectionBreakColorList = fullColorList(1:numel(sectionBreakTrials), :);


currData = inner_join(sourceData, expMd, trialMd, panelsMetadata, ftData);
currExpData = currData(strcmp(currData.expID, expID), :);
if ~isempty(trialNums)
   currExpData = currExpData(ismember(currExpData.trialNum, trialNums), :); 
end
nTrials = size(currExpData, 1);

% Put all optic flow data in a single matrix, padding with nan at ends of trials if necessary
try
flowMat = zeros(max(cellfun(@numel, currExpData.meanFlow)), nTrials);
flowFrameDur = median(diff(currExpData.frameTimes{1}));
flowFrameTimes = (1:1:size(flowMat, 1)) * flowFrameDur;
for iTrial = 1:nTrials 
    currFlow = currExpData.meanFlow{iTrial};
    currFlow(end) = median(currFlow);
    currFlow(1:20) = median(currFlow);
    if numel(currFlow) < size(flowMat, 1)
        currFlow(end + 1:size(flowMat, 1)) = nan;
    end
    flowMat(:, iTrial) = currFlow;
end
smFlowMat = repeat_smooth(flowMat, 20, 'dim', 1, 'smwin', flowSmWin); 
catch ME; rethrow(ME); end

% Put all fl data in a single matrix
flMat = zeros(currExpData.nVolumes(1), size(currExpData.rawFl{1}, 2), nTrials);
for iTrial = 1:nTrials
    flMat(:, :, iTrial) = currExpData.(flType){iTrial}; % [volume, EB wedge, trial]
end
volTimes = (1:size(flMat, 1)) ./ currExpData.volumeRate(1);

% Use optic flow data to identify movement epochs
try
moveDistSec = zeros(size(flowMat));
for iTrial = 1:nTrials
    currTrialFlow = smFlowMat(:, iTrial);
    currTrialFlow = currTrialFlow - min(currTrialFlow);
    moveFrames = currTrialFlow > moveThresh;
    moveFrameDist = zeros(1, numel(moveFrames));
    moveFrameInds = find(moveFrames);
    if ~isempty(moveFrameInds)
        for iFrame = 1:numel(moveFrames)
            moveFrameDist(iFrame) = min(abs(moveFrameInds - iFrame));
        end
    end
    moveDistSec(:, iTrial) = moveFrameDist' .* flowFrameDur';
end
catch ME; rethrow(ME); end

% Get mean panels pos data
try
meanFlShift = nan(96, size(flMat, 2), nTrials);
for iTrial = 1:nTrials
    panelsFrameTimes = double(1:currExpData.nPanelsFrames(iTrial)) ./ ...
            currExpData.panelsDisplayRate(iTrial);
    panelsPosVols = [];
    if currExpData.usingPanels(iTrial)
        for iVol = 1:size(flMat, 1)
            [~, currVol] = min(abs(panelsFrameTimes - volTimes(iVol)));
            panelsPosVols(iVol) = currExpData.panelsPosX{iTrial}(currVol);
        end
        meanFlData = [];
        flDataTemp = flMat(:, :, iTrial);
        for iPos = 1:numel(unique(currExpData.panelsPosX{iTrial}))
            meanFlData(iPos, :) = ...
                mean(flDataTemp(panelsPosVols == (iPos - 1), :), 1, 'omitnan'); % --> [barPos, wedge]
        end
        meanFlShift(:, :, iTrial) = cat(1, meanFlData(92:96, :), meanFlData(1:91, :));
        cm = hsv(size(meanFlShift, 2)) .* 0.9;
    end
end
catch ME; rethrow(ME); end

disp('Block data ready')

%% PLOT SUMMARY OF MOVEMENT THROUGHOUT EXPERIMENT

try 
    
% Heatmap of slightly smoothed flow throughout experiment
f = figure(1);clf;
f.Color = [1 1 1];
ax = gca;


imagesc([0, currExpData.trialDuration(1)], [1, nTrials], smFlowMat')
colormap(viridis)
title('Optic flow')
ax.YTick = 1:nTrials;
xlabel('Time (sec)')
ylabel('Trial')
ax.FontSize = 12;

% Plot avg flow value and percentage of movement vols for each trial
f = figure(2); clf; hold on;
f.Color = [1 1 1];
f.Position(3:4) = [1200 400];
ax = subaxis(1, 1, 1, 'mt', 0.1, 'mb', 0.16, 'ml', 0.08, 'mr', 0.08); hold on
avgTrialFlow = mean(smFlowMat, 1, 'omitnan');
plotX = 1:numel(avgTrialFlow);
plot(plotX, avgTrialFlow, '-o', 'color', 'b', 'linewidth', 1, 'markersize', 8)
ax.YColor = [0 0 1];
ax.FontSize = 12;
xlabel('Trial', 'fontsize', 16)
xlim([0 plotX(end) + 1])
ax.XTick = 1:nTrials;
ax.XTickLabel = currExpData.trialNum;

% Plot lines dividing sections if necessary
yL = ylim();
if ~isempty(sectionBreakTrials)
    for iColor = 1:size(sectionBreakColorList, 1)
        currColor = sectionBreakColorList(iColor, :);
        currSectionBreaks = sectionBreakTrials(sectionBreakColors(iColor));
        for iSection = 1:numel(currSectionBreaks)
            currTrialNum = currSectionBreaks(iSection);
            xx = [currTrialNum currTrialNum] - 0.5;
            yy = yL + [-100 100];
            if iColor == 1 && iSection == 1
                plotLines = plot(xx, yy, 'color', currColor, 'linewidth', 2);
            elseif iSection == 1
                plotLines(end + 1) = plot(xx, yy, 'color', currColor, 'linewidth', 2);
            else
                plot(xx, yy, 'color', currColor, 'linewidth', 2);
            end
        end
    end
    legend(plotLines, sectionColorLabels, 'autoupdate', 'off', 'location', 'best')
end
ylim(yL)

% Plot proportion of volumes when fly was moving in each trial
plotData = smFlowMat - min(smFlowMat(:));
moveFrames = smFlowMat > moveThresh;
trialMovePercent = (sum(moveFrames) ./ size(moveFrames, 1)) * 100;
ylabel('Mean optic flow (AU)', 'fontsize', 16)
yyaxis('right');
ax.YColor = [1 0 0];
ylabel('% movement', 'fontsize', 16)
plot(trialMovePercent, '-*', 'color', 'r', 'linewidth', 1, 'markersize', 10)
ylim([0 100])
yyaxis('left');
ax.YColor = [0 0 1];

catch ME; rethrow(ME); end

%% PLOT SINGLE-TRIAL HEATMAPS FOR ENTIRE BLOCK

smWin = 5;

omitMoveVols = 0;
moveDistThresh = 2;

figPos = [];

try 
    
plotFl = smoothdata(flMat, 1, 'gaussian', smWin, 'omitnan');

% Replace volumes when fly was moving with nan if necessary
if omitMoveVols
   moveDistArr = permute(repmat(moveDistSec, 1, 1, size(flMat, 2)), [1 3 2]);
   plotFl(moveDistArr < moveDistThresh) = nan; 
end

% To give all figures the same color scale
plotFl(1, 1, :) = min(plotFl(:), [], 'omitnan');
plotFl(end, end, :) = max(plotFl(:), [], 'omitnan');

f = figure(10);clf;
f.Color = [1 1 1];
if ~isempty(figPos)
    f.Position(3:4) = figPos;
end
colormap(magma);
for iTrial = 1:nTrials 
    
    subaxis(nTrials, 8, [1:7] + (8 * (iTrial - 1)), 'mt', 0.01, 'mb', 0.11, 'sv', 0.002, 'mr', 0.02, ...
            'ml', 0.04, 'sh', 0.02);
    imagesc([0, currExpData.trialDuration(1)], [1, 8], plotFl(:, :, iTrial)');
    hold on
    plot(volTimes, 9 - smoothdata(currExpData.pvaWedge{iTrial}, 1, 'gaussian', smWin), 'color', 'g')
    
    % Label X axis on final trial
    ax = gca();
    if iTrial < nTrials
        ax.XTickLabel = [];
    else
        ax.FontSize = 12;
        xlabel('Time (sec)', 'FontSize', 14);
    end
    
    % Label each plot with trial number
    ax.YTickLabel = [];
    t = ylabel(['Trial #', num2str(currExpData.trialNum(iTrial))], 'rotation', 0, 'FontSize', 12);
    t.HorizontalAlignment = 'right';
    t.VerticalAlignment = 'middle';
    t.Position(1) = t.Position(1) * 2;
    
    % Add fly movement summary plot to the side of each trial
    ax = subaxis(nTrials, 8, 8*iTrial);
    b = barh(2, trialMovePercent(iTrial) / 100, 'facecolor', [0.75 0.1 0.1]);
    ax.YLim = [0.5 3.5];
    ax.XLim = [0 1];
    ax.YTickLabel = [];
    
    if iTrial < nTrials
        ax.XTickLabel = [];
    else
        ax.XTick = [0 1];
        ax.XTickLabel = [0 100];
        xlabel('% movement')
        ax.FontSize = 12;
    end
    
end

catch ME; rethrow(ME); end

%% PLOT TUNING HEATMAPS FOR EACH WEDGE ACROSS TRIALS

smWin = 5;

omitMoveVols = 0;
moveDistThresh = 2;

figPos = [];

try 
    
plotFl = smoothdata(meanFlShift, 1, 'gaussian', smWin, 'omitnan');

% Determine subplot grid size
nPlots = size(plotFl, 2);
if nPlots == 3
    plotPos = [2 2]; % Because [1 3] looks bad
else
    plotPos = numSubplots(nPlots);
end

% Create figure and plots
f = figure(1);clf;
f.Color = [1 1 1];
for iPlot = 1:nPlots
    
    subaxis(plotPos(1), plotPos(2), iPlot, 'ml', 0.04, 'mr', 0.03, 'mb', 0.1);
    
    currData = squeeze(plotFl(:, iPlot, :)); % --> [barPos, trial]
    currDataShift = [currData(92:96, :); currData(1:91, :)];
    
    % Plot data
    plotX = -180:3.75:(180 - 3.75);
    imagesc(plotX, [1 size(currData, 2)], currDataShift');
    hold on; plot([0 0], [0 size(currData, 2) + 1], 'color', 'k', 'linewidth', 2)
    ylim([0.5 size(currData, 2) + 0.5])
    xlabel('Bar position (degrees from front of fly)');
    ylabel('Trial');
    colormap('magma')
    
    % Add Y tick labels
    ax = gca;
    ax.YTick = 1:size(currData, 2);
    ax.YTickLabels = currExpData.trialNum;
    
    ax.FontSize = 10
    
end%iPlot

catch ME; rethrow(ME); end

%% PLOT AS LINES INSTEAD OF USING IMAGESC

smWin = 5;

omitMoveVols = 0;
moveDistThresh = 2;

figPos = [];

try 
    
plotFl = smoothdata(meanFlShift, 1, 'gaussian', smWin, 'omitnan');

% Offset data so that each plot is centered at zero
shiftDataOffset = [];
for iTrial = 1:size(plotFl, 3)
   for iWedge = 1:size(plotFl, 2)
      currData = plotFl(:, iWedge, iTrial); % --> [barPos]
      shiftDataOffset(:, iWedge, iTrial) = currData - mean(currData, 'omitnan');
   end    
end

% Find max and min values 
yMax = max(shiftDataOffset(:), [], 'omitnan');
yMin = min(shiftDataOffset(:), [], 'omitnan');
range = max(abs([yMax, yMin]), [], 'omitnan') *  2.5;

% Determine subplot grid size
nPlots = size(plotFl, 2);
if nPlots == 3
    plotPos = [2 2]; % Because [1 3] looks bad
else
    plotPos = numSubplots(nPlots);
end

% Create figure and plots
f = figure(1);clf;
f.Color = [1 1 1];

allAx = {}; plotLines = [];
for iPlot = 1:size(shiftDataOffset, 2)
    allAx{iPlot} = subaxis(1, nPlots, iPlot, 'mt', 0, 'ml', 0.05, 'mr', 0.03, 'mb', 0.08, 'sh', 0.05);
    hold on;
    currData = squeeze(shiftDataOffset(:, iPlot, :)); % --> [barPos, trial]

    % Plot data
    plotX = -180:3.75:(180 - 3.75);
    
    % Separate data into distinct rows
    offsets = []; allPlotX = [];
    currData = fliplr(currData); % Flip so first trial is in the last column
    for iTrial = 1:size(currData, 2)
        currOffset = (range * (iTrial - 1));
        offsets(iTrial) = currOffset;
        currData(:, iTrial) = currData(:, iTrial) + currOffset;
        allPlotX(:, iTrial) = plotX;
        allPlotX(isnan(currData(:, iTrial)), iTrial) = nan;
%         plot([plotX(1), plotX(end)], [currOffset, currOffset], '--', 'color', 'b') % Plot zero lines
    end
    plot(allPlotX, currData, 'color', 'k', 'linewidth', 1.25);
    yL(1) = yMin * 1.25;
    yL(2) = range * (size(currData, 2));
    hold on; plot([0 0], yL, 'color', 'b', 'linewidth', 1)
    ylim(yL);
    
    xlabel('Bar position (deg)');
    if iPlot == 1
        ylabel('Trial', 'fontsize', 16);
    end
    allAx{iPlot}.YTick = offsets;
    yTickLabel = size(shiftDataOffset, 3):-1:1;
    
end%iPlot

catch ME; rethrow(ME); end

