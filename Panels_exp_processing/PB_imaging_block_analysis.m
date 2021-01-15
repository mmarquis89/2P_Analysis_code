
%% Load data for all experiments
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData_60D05_7f';
figDir = fullfile(parentDir, 'Figs');

expList = [{'20201114-1', '20201117-1', '20201117-2', '20201120-2', '20201201-1', '20201201-2', ...
        '20201201-3', '20201203-1', '20201203-2', '20201210-1'}];%%

[expMd, trialMd, roiData, ftData, flailingEvents, panelsMetadata, wedgeData, glomData] = ...
        load_PB_data(parentDir, expList);

% Load file with info about the details of the bath applications
drugTimingMd = readtable(fullfile(parentDir, 'drug_timing_data.csv'), 'delimiter', ',');

%% Group data from several compatible trials into a single block

% expID = '20201210-1';
expID = expMd.expID{9};

trialNums = [];

sourceData = wedgeData; %glomData;%

showBehaviorPlots = 0;

washoutTime = 0;

flSmWin = 5;

flowSmWin = 6;
moveThresh = 0.08;%unique(flailingEvents.eventData.moveThresh(strcmp(flailingEvents.eventData.expID, ...
        %expID)));

try 
    
% Generate consolidated source data table
currData = inner_join(sourceData, expMd, trialMd, panelsMetadata);
currData = outerjoin(currData, ftData, 'type', 'left', 'mergekeys', 1);
currExpData = currData(strcmp(currData.expID, expID), :);
currExpData = outerjoin(currExpData, drugTimingMd, 'type', 'left', 'mergekeys', 1);
if ~isempty(trialNums)
   currExpData = currExpData(ismember(currExpData.trialNum, trialNums), :); 
end
nTrials = size(currExpData, 1);
nRois = size(currExpData.rawFl{1}, 2);

% Put all optic flow data in a single matrix, padding with nan at ends of trials if necessary
try
flowMat = zeros(max(cellfun(@numel, currExpData.meanFlow)), nTrials);
for iTrial = 1:nTrials 
    currFlow = currExpData.meanFlow{iTrial};
    if ~isempty(currFlow)
        currFlow(end) = median(currFlow);
        currFlow(1:20) = median(currFlow);
        if numel(currFlow) < size(flowMat, 1)
            currFlow(end + 1:size(flowMat, 1)) = nan;
        end
        flowMat(:, iTrial) = currFlow;
    end
end
smFlowMat = repeat_smooth(flowMat, 20, 'dim', 1, 'smwin', flowSmWin); 
catch ME; rethrow(ME); end

% Put all fl data in a single matrix for each fl type
rawFlMat = zeros(currExpData.nVolumes(1), nRois, nTrials);
trialDffMat = rawFlMat;
expDffMat = rawFlMat;
for iTrial = 1:nTrials
    rawFlMat(:, :, iTrial) = currExpData.rawFl{iTrial}; % [volume, ROI, trial]
    trialDffMat(:, :, iTrial) = currExpData.trialDff{iTrial}; % [volume, ROI, trial]
    expDffMat(:, :, iTrial) = currExpData.expDff{iTrial}; % [volume, ROI, trial]
end
volTimes = (1:size(rawFlMat, 1)) ./ currExpData.volumeRate(1);
nVolumes = numel(volTimes);

% Use optic flow data to identify movement epochs
try
% moveDistSec = zeros(size(flowMat));
for iTrial = 1:nTrials
    currTrialFlow = smFlowMat(:, iTrial);
    currTrialFlow = currTrialFlow - min(currTrialFlow);
    moveFrames = currTrialFlow > moveThresh;
    moveFrameDist = zeros(1, numel(moveFrames));
    moveFrameInds = find(moveFrames);
%     if ~isempty(moveFrameInds)
%         for iFrame = 1:numel(moveFrames)
%             moveFrameDist(iFrame) = min(abs(moveFrameInds - iFrame));
%         end
%     end
%     moveDistSec(:, iTrial) = moveFrameDist' .* flowFrameDur';
end
catch ME; rethrow(ME); end

% Get mean fluorescence at each bar location for all ROIs
try
rawFlTuningShifted = nan(96, nRois, nTrials);
trialDffTuningShifted = nan(96, nRois, nTrials);
expDffTuningShifted = nan(96, nRois, nTrials);
drugTrials = currExpData.trialNum(~isnan(currExpData.startTime));
for iTrial = 1:nTrials
    panelsFrameTimes = double(1:currExpData.nPanelsFrames(iTrial)) ./ ...
            currExpData.panelsDisplayRate(iTrial);
    panelsPosVols = [];
    if currExpData.usingPanels(iTrial)
        
        % Identify bar location during each volume
        for iVol = 1:nVolumes
            [~, currVol] = min(abs(panelsFrameTimes - volTimes(iVol)));
            panelsPosVols(iVol) = currExpData.panelsPosX{iTrial}(currVol);
        end
        
        % Get Fl data for current trial
        rawFlTuningData = [];
        trialDffTuningData = [];
        expDffTuningData = [];
        currRawFl = rawFlMat(:, :, iTrial);
        currTrialDff = trialDffMat(:, :, iTrial);
        currExpDff = expDffMat(:, :, iTrial);
        
        % Set volumes before and during drug application to nan
        if ismember(currExpData.trialNum(iTrial), drugTrials)
            currRawFl(volTimes < (currExpData.startTime(iTrial) + ...
                    currExpData.duration(iTrial) + washoutTime), :) = nan;
            currTrialDff(volTimes < (currExpData.startTime(iTrial) + ...
                    currExpData.duration(iTrial) + washoutTime), :) = nan;
            currExpDff(volTimes < (currExpData.startTime(iTrial) + ...
                    currExpData.duration(iTrial) + washoutTime), :) = nan;
        end
        
        % Calculate mean bar position tuning
        for iPos = 1:numel(unique(currExpData.panelsPosX{iTrial}))
            rawFlTuningData(iPos, :) = ...
                mean(currRawFl(panelsPosVols == (iPos - 1), :), 1, 'omitnan'); % --> [barPos, ROI]
            trialDffTuningData(iPos, :) = ...
                mean(currTrialDff(panelsPosVols == (iPos - 1), :), 1, 'omitnan'); % --> [barPos, ROI]
            expDffTuningData(iPos, :) = ...
                mean(currExpDff(panelsPosVols == (iPos - 1), :), 1, 'omitnan'); % --> [barPos, ROI]
        end
        
        % Shift data so it is centered around the bar at 0°
        rawFlTuningShifted(:, :, iTrial) = cat(1, rawFlTuningData(92:96, :), ...
                rawFlTuningData(1:91, :));
        trialDffTuningShifted(:, :, iTrial) = cat(1, trialDffTuningData(92:96, :), ...
                trialDffTuningData(1:91, :));
        expDffTuningShifted(:, :, iTrial) = cat(1, expDffTuningData(92:96, :), ...
                expDffTuningData(1:91, :));
        cm = hsv(nRois) .* 0.9;
    end
end

catch ME; rethrow(ME); end

% Calculate whole-trial vector strength and phase (skipping drug applications)
try
fullTrialVectorStrength = [];
fullTrialVectorPhase = [];
for iTrial = 1:nTrials
    
    if currExpData.usingPanels(iTrial)
        panelsFrameTimes = double(1:currExpData.nPanelsFrames(iTrial)) ./ ...
                currExpData.panelsDisplayRate(iTrial);
        panelsPosVols = [];
        
        % Identify bar location during each volume
        for iVol = 1:nVolumes
            [~, currVol] = min(abs(panelsFrameTimes - volTimes(iVol)));
            panelsPosVols(iVol) = currExpData.panelsPosX{iTrial}(currVol);
        end
        
        % Get Fl data for current trial
        currExpDff = expDffMat(:, :, iTrial);
        
        % Set volumes before and during drug application to nan
        if ismember(currExpData.trialNum(iTrial), drugTrials)
            currExpDff(volTimes < (currExpData.startTime(iTrial) + ...
                    currExpData.duration(iTrial) + washoutTime), :) = nan;
        end
        
        % Convert panels position to a phase
        panelsPhase = 2 * pi * (panelsPosVols ./ max(panelsPosVols));
        
        % Calculate expDff vector strength and phase for current trial
        for iRoi = 1:nRois
            currFlData = smoothdata(currExpDff(:, iRoi), 1, 'gaussian', flSmWin);
            currFlData = currFlData(~isnan(currFlData));
            currPanelsPhase = panelsPhase(~isnan(currFlData));
            oppositeVector = sum(currFlData' .* sin(currPanelsPhase));
            adjacentVector = sum(currFlData' .* cos(currPanelsPhase));
            oppositeMeanVector =  oppositeVector / numel(currPanelsPhase);
            adjacentMeanVector =  adjacentVector / numel(currPanelsPhase);
            
            fullTrialVectorStrength(iTrial, iRoi) = sqrt(oppositeVector^2 + adjacentVector^2) / ...
                    numel(currPanelsPhase);
            fullTrialVectorPhase(iTrial, iRoi) = wrapTo2Pi(atan(oppositeMeanVector / ...
                    adjacentMeanVector));
            
        end
    end
end
catch ME; rethrow(ME); end

% Calculate vector strength/phase and extract dF/F data for each individual bar rotation cycle
try
cycleVectorStrength = [];
cycleVectorPhase = [];
cycleFlData = [];
cycleDffData = {};
allCycleStartVols = {};
cycleVolTimes = [];
for iTrial = 1:nTrials
    
    panelsFrameTimes = double(1:currExpData.nPanelsFrames(iTrial)) ./ ...
        currExpData.panelsDisplayRate(iTrial);
    if currExpData.usingPanels(iTrial)
        panelsPosVols = [];
        
        % Identify bar location during each volume
        for iVol = 1:nVolumes
            [~, currVol] = min(abs(panelsFrameTimes - volTimes(iVol)));
            panelsPosVols(iVol) = currExpData.panelsPosX{iTrial}(currVol);
        end
        
        % Identify the start and end of each bar cycle
        cycleStartVolsStr = regexprep(num2str(panelsPosVols == 0), ' ', '');
        cycleStartVols = regexp(cycleStartVolsStr, '01');
        cycleStartVols = [1, cycleStartVols];
        cycleEndVols = cycleStartVols(2:end) - 1;
        cycleEndVols = [cycleEndVols, nVolumes];
        nCycles = numel(cycleStartVols);
        allCycleStartVols{iTrial} = cycleStartVols;
        
        % Get Fl data for current trial
        currExpDff = expDffMat(:, :, iTrial);
        
        % Convert panels position to a phase
        panelsPhase = 2 * pi * (panelsPosVols ./ max(panelsPosVols));
        
        % Calculate vector strength and phase
        trialCycleVolTimes = [];
        for iRoi = 1:nRois
            for iCycle = 1:nCycles
                startInd = cycleStartVols(iCycle);
                endInd = cycleEndVols(iCycle);
                currFlData = smoothdata(currExpDff(startInd:endInd, iRoi), 1, 'gaussian', flSmWin);
                currVolTimes = volTimes(cycleStartVols(iCycle):cycleEndVols(iCycle))';
                currPanelsPhase = panelsPhase(startInd:endInd);
                cycleFlData{iTrial, iRoi, iCycle} = currFlData;
                trialCycleVolTimes{iCycle} = currVolTimes;
                oppositeVector = sum(currFlData' .* sin(currPanelsPhase));
                adjacentVector = sum(currFlData' .* cos(currPanelsPhase));
                oppositeMeanVector =  oppositeVector / numel(currPanelsPhase);
                adjacentMeanVector =  adjacentVector / numel(currPanelsPhase);
                
                cycleVectorStrength(iTrial, iRoi, iCycle) = sqrt(oppositeVector^2 ...
                        + adjacentVector^2) / numel(currPanelsPhase);
                cycleVectorPhase(iTrial, iRoi, iCycle) = wrapTo2Pi(atan2(oppositeMeanVector, ...
                        adjacentMeanVector));
            end
        end
        
        % Extract dF/F data, pad with nan, and consolidate into single array
        currTrialCycleFl = squeeze(cycleFlData(iTrial, :, :));
        maxVolCount = max(cellfun(@numel, currTrialCycleFl(1, :)));
        currCycleDffArr = [];
        currCycleVolTimesArr = [];
        for iCycle = 1:nCycles
            currCycleFl = cell2mat(currTrialCycleFl(:, iCycle)');
            currCycleVolTimes = trialCycleVolTimes{iCycle};
            if size(currCycleFl, 1) < maxVolCount
                currCycleFl(size(currCycleFl, 1) + 1:maxVolCount, :) = nan;
            end
            if numel(currCycleVolTimes) < maxVolCount
                currCycleVolTimes(numel(currCycleVolTimes) + 1:maxVolCount) = nan;
            end
            currCycleDffArr = cat(3, currCycleDffArr, currCycleFl); % [vol, glom, cycle]
            currCycleVolTimesArr = cat(2, currCycleVolTimesArr, currCycleVolTimes);
        end
        cycleDffData{iTrial} = currCycleDffArr;
        cycleVolTimes{iTrial} = currCycleVolTimesArr;
    end
end
catch ME; rethrow(ME); end

disp('Block data ready')

catch ME; rethrow(ME); end

% PLOT SUMMARY OF MOVEMENT THROUGHOUT EXPERIMENT

try 
% Heatmap of slightly smoothed flow throughout experiment
f = figure(1);clf;
f.Color = [1 1 1];
ax = gca;
imagesc([0, currExpData.trialDuration(1)], [1, nTrials], smFlowMat')
colormap(viridis)
hold on;
drugTrials = find(~isnan(currExpData.startTime));
if ~isempty(drugTrials)
    for iTrial = 1:numel(drugTrials)
        if strcmp(currExpData.visualStim{drugTrials(iTrial)}, 'yes')
            lineColor = 'r';
        else
            lineColor = 'm';
        end
        xx_1 = [1, 1] .* currExpData.startTime(drugTrials(iTrial));
        xx_2 = [1, 1] .* currExpData.startTime(drugTrials(iTrial)) + ...
                currExpData.duration(drugTrials(iTrial));
        yy = [-0.5, 0.5] + currExpData.trialNum(drugTrials(iTrial));
        plot(xx_1, yy, 'color', lineColor, 'linewidth', 5);
        plot(xx_2, yy, 'color', lineColor, 'linewidth', 5);
    end
end
titleStr = {[expID, '  —  Optic flow']};
if ~isempty(drugTrials)
    titleStr = [titleStr, {['(red lines = ', currExpData.drugName{drugTrials(iTrial)}, ...
            ' onset and offset)']}];
end
title(titleStr);
ax.YTick = 1:nTrials;
xlabel('Time (sec)')
ylabel('Trial')
ax.FontSize = 14;
if ~showBehaviorPlots
   close(f); 
end

% Plot avg flow value and percentage of movement vols for each trial (after the onset of any drug
% application)
f = figure(2); clf; hold on;
f.Color = [1 1 1];
f.Position(3:4) = [1000 350];
ax = subaxis(1, 1, 1, 'mt', 0.1, 'mb', 0.16, 'ml', 0.08, 'mr', 0.08); hold on
postDrugFlowMat = smFlowMat;
if ~isempty(drugTrials)
    for iTrial = 1:numel(drugTrials)
        currTrialNum = drugTrials(iTrial);
        postDrugFlowMat(currExpData.vidFrameTimes{currTrialNum} < ...
                currExpData.startTime(currTrialNum)) = nan;
    end
end
avgTrialFlow = mean(postDrugFlowMat, 1, 'omitnan');
xx = 1:numel(avgTrialFlow);
plot(xx, avgTrialFlow, '-o', 'color', 'b', 'linewidth', 1, 'markersize', 8)
ax.YColor = [0 0 1];
ax.FontSize = 12;
xlabel('Trial', 'fontsize', 16)
xlim([0 xx(end) + 1])
ax.XTick = 1:nTrials;
ax.XTickLabel = currExpData.trialNum;
ylabel('Mean optic flow (AU)', 'fontsize', 16)
ylim(ylim() .* [0.95, 1.05]);

% Plot proportion of volumes when fly was moving in each trial
plotData = postDrugFlowMat - min(postDrugFlowMat(:), [], 'omitnan');
moveFrames = postDrugFlowMat > moveThresh;
trialMovePercent = [];
for iTrial = 1:size(moveFrames, 2)
    currMoveFrames = moveFrames(:, iTrial);
    trialMovePercent(iTrial) = (sum(currMoveFrames, 'omitnan') ./ ...
            sum(~isnan(currMoveFrames))) * 100;
end
yyaxis('right');
ax.YColor = [1 0 0];
ylabel('% movement', 'fontsize', 16)
plot(trialMovePercent, '-*', 'color', 'r', 'linewidth', 1, 'markersize', 10)
hold on;
if ~isempty(drugTrials)
    plot_stim_shading([drugTrials - 0.25, drugTrials + 0.25], 'color', [0 0.9 0])
end
ylim([0 100])
yyaxis('left');
ax.YColor = [0 0 1];
titleStr = [expID, '  —  Fly movement summary'];
if ~isempty(drugTrials)
    titleStr = [titleStr, ' (green = ', currExpData.drugName{drugTrials(1)}, ' trials)'];  
end
t = title(titleStr);
t.Position(2) = t.Position(2) * 1.01;
xlim(xlim() + [0.5, -0.5])
box off
if ~showBehaviorPlots
   close(f)
end

catch ME; rethrow(ME); end

%% PLOT SINGLE-TRIAL HEATMAPS FOR ENTIRE BLOCK

saveFig = 0;

smWin = 5;
flType = 'expDff';
plotBarCycles = 1;

figPos = [1800 950];

% Adjust plot spacing and margins
SV = 0.002;
SH = 0.02;
ML = 0.04;
MR = 0.02;
MT = 0.05;
MB = 0.065;

try 
    
% Get correct source data
if strcmp(flType, 'rawFl')
    flData = rawFlMat;
elseif strcmp(flType, 'trialDff')
    flData = trialDffMat;
elseif strcmp(flType, 'expDff')
    flData = expDffMat;
end
plotFl = smoothdata(flData, 1, 'gaussian', smWin, 'omitnan');

% To give all figures the same color scale
plotFl(1, 1, :) = min(plotFl(:), [], 'omitnan');
plotFl(end, end, :) = max(plotFl(:), [], 'omitnan');

f = figure(3);clf;
f.Color = [1 1 1];
if ~isempty(figPos)
    f.Position(3:4) = figPos;
    if sum(f.Position([2, 4])) > 1080
        f.Position(2) = 50;
    end
end
colormap(magma);
drugTrials = find(~isnan(currExpData.startTime));
for iTrial = 1:nTrials 
    
    % Plot Fl Data and PVA
    subaxis(nTrials, 8, [1:7] + (8 * (iTrial - 1)), 'mt', MT, 'mb', MB, 'sv', SV, 'mr', MR, ...
            'ml', ML, 'sh', SH);
    imagesc([0, currExpData.trialDuration(1)], [1, 8], plotFl(:, :, iTrial)');
    hold on
    if size(plotFl, 2) == 8
        plot(volTimes, 9 - smoothdata(currExpData.pvaWedge{iTrial}, 1, 'gaussian', smWin), ...
                'color', 'c', 'linewidth', 0.5)
    end
    
    % Plot any drug applications in the current trial
    if ~isempty(drugTrials) && ismember(iTrial, drugTrials)
        if strcmp(currExpData.visualStim{iTrial}, 'yes')
            lineColor = 'g';
        else
            lineColor = rgb('yellow');
        end
        yL = ylim();
        xx_1 = [1, 1] .* currExpData.startTime(iTrial);
        xx_2 = [1, 1] .* currExpData.startTime(iTrial) + currExpData.duration(iTrial);
        plot(xx_1, yL, 'color', lineColor, 'linewidth', 5);
        plot(xx_2, yL, 'color', lineColor, 'linewidth', 5);
        ylim(yL)
    end
    
    % Plot bar cycle boundaries if necessary
    currStartVols = allCycleStartVols{iTrial};
    if plotBarCycles && ~isempty(currStartVols)
        yL = ylim();
        for iCycle = 1:numel(currStartVols)
            xx = [1, 1] .* volTimes(currStartVols(iCycle));
            plot(xx, yL, 'color', 'r', 'linewidth', 1)
        end
        ylim(yL);
    end
    
    % Label each plot with trial number
    ax = gca();
    ax.YTickLabel = [];
    t = ylabel(['Trial #', num2str(currExpData.trialNum(iTrial))], 'rotation', 0, 'FontSize', 12);
    t.HorizontalAlignment = 'right';
    t.VerticalAlignment = 'middle';
    t.Position(1) = t.Position(1) * 2;
    
    % Label X axis on final trial
    if iTrial < nTrials
        ax.XTickLabel = [];
    else
        str = ax.YLabel.String;
        pos = ax.YLabel.Position;
        ax.FontSize = 12;
        xlabel('Time (sec)', 'FontSize', 14);
        ylabel(str, 'rotation', 0, 'FontSize', 12);      
    end
    
    % Add title above the first plot
    if iTrial == 1
        titleStr = [expID, '  —  single trial ', flType];
        if ~isempty(drugTrials)
            titleStr = [titleStr, '  (green lines = ', currExpData.drugName{drugTrials(1)}, ...
                ' onset and offset)'];
        end
        title(titleStr, 'FontSize', 16);
    end
    
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

% Save figure
if saveFig
    figTitle = [expID, '_single_trial_heatmaps'];
    save_figure(f, figDir, figTitle);
end

catch ME; rethrow(ME); end

%% PLOT TUNING HEATMAPS FOR EACH WEDGE ACROSS TRIALS

saveFig = 0;

smWin = 5;

flType = 'expDff';

figPos = [1600 950];

% Adjust plot spacing and margins
SV = 0.01;
SH = 0.03;
ML = 0.04;
MR = 0.02;
MT = 0;
MB = 0.08;

try 

% Get correct source data
if strcmp(flType, 'rawFl')
    flData = rawFlTuningShifted;
elseif strcmp(flType, 'trialDff')
    flData = trialDffTuningShifted;
elseif strcmp(flType, 'expDff')
    flData = expDffTuningShifted;
end
plotFl = smoothdata(flData, 1, 'gaussian', smWin, 'omitnan');

% Determine subplot grid size
nPlots = size(plotFl, 2);
if nPlots == 3
    plotPos = [2 2]; % Because [1 3] looks bad
else
    plotPos = numSubplots(nPlots);
end

% Create figure and plots
f = figure(4);clf;
if ~isempty(figPos)
    f.Position(3:4) = figPos;
    if sum(f.Position([2, 4])) > 1080
        f.Position(2) = 50;
    end
end
f.Color = [1 1 1];
for iPlot = 1:nPlots
    
    subaxis(plotPos(1), plotPos(2), iPlot, 'mt', MT, 'mb', MB, 'sv', SV, 'mr', MR, ...
            'ml', ML, 'sh', SH);
    
    currData = squeeze(plotFl(:, iPlot, :)); % --> [barPos, trial]
    
    % Plot data
    xx = -180:3.75:(180 - 3.75);
    imagesc(xx, [1 size(currData, 2)], currData');
    hold on; plot([0 0], [0 size(currData, 2) + 1], 'color', 'k', 'linewidth', 2)
    ylim([0.5 size(currData, 2) + 0.5])
    xlabel('Bar position (degrees from front of fly)');
    ylabel('Trial');
    colormap('magma')
    
    % Add Y tick labels
    ax = gca;
    ax.YTick = 1:size(currData, 2);
    ax.YTickLabels = currExpData.trialNum;
    
    % Remove extraneous X and Y axis labels
    if mod(iPlot - 1, plotPos(2))
        ax.YLabel.String = [];
    end
    if iPlot <= (nPlots - plotPos(2))
        ax.XLabel.String = [];
        ax.XTickLabel = [];
    end
    
    ax.FontSize = 12;
    
    % Indicate transitions to drug application trials
    if ~isempty(drugTrials)
        for iTrial = 1:numel(drugTrials)
            xL = xlim();
            if strcmp(currExpData.visualStim{drugTrials(iTrial)}, 'yes')
                lineColor = 'g';
            else
                lineColor = rgb('yellow');
            end
            plot(xL, [-0.5, -0.5] + drugTrials(iTrial), 'color', lineColor, 'linewidth', 3)
        end
    end
    
end%iPlot

% Add title to figure
titleStr = [expID, '  —  EB wedge mean ', flType, ' visual tuning'];
if ~isempty(drugTrials)
    titleStr = [titleStr, '  (green lines = ', currExpData.drugName{drugTrials(1)}, ' application)'];
end
h = suptitle(titleStr);
h.FontSize = 18;
    
% Save figure
if saveFig
    figTitle = [expID, '_visual_tuning_heatmaps'];
    save_figure(f, figDir, figTitle);
end

catch ME; rethrow(ME); end


%% PLOT AS LINES INSTEAD OF USING IMAGESC

saveFig = 0;

smWin = 5;

flType = 'expDff';

figPos = [1700 950];

% Adjust plot spacing and margins
SV = 0.01;
SH = 0.03;
ML = 0.04;
MR = 0.02;
MT = 0;
MB = 0.08;

try 
    
% Get correct source data
if strcmp(flType, 'rawFl')
    flData = rawFlTuningShifted;
elseif strcmp(flType, 'trialDff')
    flData = trialDffTuningShifted;
elseif strcmp(flType, 'expDff')
    flData = expDffTuningShifted;
end
plotFl = smoothdata(flData, 1, 'gaussian', smWin, 'omitnan');

% Offset data so that each plot is centered at zero
flDataOffset = [];
for iTrial = 1:size(plotFl, 3)
   for iWedge = 1:size(plotFl, 2)
      currData = plotFl(:, iWedge, iTrial); % --> [barPos]
      flDataOffset(:, iWedge, iTrial) = currData - mean(currData, 'omitnan');
   end    
end

% Find max and min values 
yMax = max(flDataOffset(:), [], 'omitnan');
yMin = min(flDataOffset(:), [], 'omitnan');
range = max(abs([yMax, yMin]), [], 'omitnan') *  2.5;

% Determine subplot grid size
nPlots = size(plotFl, 2);
if nPlots == 3
    plotPos = [2 2]; % Because [1 3] looks bad
else
    plotPos = numSubplots(nPlots);
end

% Create figure and plots
f = figure(5);clf;
f.Color = [1 1 1];
if ~isempty(figPos)
    f.Position(3:4) = figPos;
    if sum(f.Position([2, 4])) > 1080
        f.Position(2) = 50;
    end
end
allAx = {}; plotLines = [];
for iPlot = 1:size(flDataOffset, 2)
    allAx{iPlot} = subaxis(1, nPlots, iPlot, 'mt', MT, 'mb', MB, 'sv', SV, 'mr', MR, ...
            'ml', ML, 'sh', SH);
    hold on;
    currData = squeeze(flDataOffset(:, iPlot, :)); % --> [barPos, trial]

    % Plot data
    xx = -180:3.75:(180 - 3.75);
    
    % Separate data into distinct rows
    offsets = []; allPlotX = [];
    currData = fliplr(currData); % Flip so first trial is in the last column
    for iTrial = 1:size(currData, 2)
        currOffset = (range * (iTrial - 1));
        offsets(iTrial) = currOffset;
        currData(:, iTrial) = currData(:, iTrial) + currOffset;
        allPlotX(:, iTrial) = xx;
        allPlotX(isnan(currData(:, iTrial)), iTrial) = nan;
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
    allAx{iPlot}.YTickLabel = size(flDataOffset, 3):-1:1;
    allAx{iPlot}.FontSize = 14;
    
    % Indicate transitions to drug application trials
    if ~isempty(drugTrials)
        offsetsRev = offsets(end:-1:1);
        for iTrial = 1:numel(drugTrials)
            xL = xlim();
            if strcmp(currExpData.visualStim{drugTrials(iTrial)}, 'yes')
                lineColor = 'g';
            else
                lineColor = rgb('yellow');
            end
            plotOffset = offsetsRev(drugTrials(iTrial) - 1) - (range/2);
            plot(xL, [1, 1] .* plotOffset, 'color', lineColor, 'linewidth', 3)
        end
    end
    
end%iPlot

titleStr = [expID, '  —  EB wedge mean ', flType, ' visual tuning'];
if ~isempty(drugTrials)
    titleStr = [titleStr, '  (green lines = ', currExpData.drugName{drugTrials(1)}, ' application)'];
end
h = suptitle(titleStr);
h.FontSize = 18;

% Save figure
if saveFig
    figTitle = [expID, '_visual_tuning_curves'];
    save_figure(f, figDir, figTitle);
end

catch ME; rethrow(ME); end


%% PLOT MIN AND MAX VALUES FROM THE TUNING CURVES

saveFig = 0;

smWin = 5;

flType = 'expDff';
plotType = 'minmax'; % minmax, amplitude

figPos = [550 950];

% Adjust plot spacing and margins
SV = 0.03;
SH = 0.03;
ML = 0.1;
MR = 0.02;
MT = 0.1;
MB = 0.07;

try

% Get correct source data
if strcmp(flType, 'rawFl')
    flData = rawFlTuningShifted;
elseif strcmp(flType, 'trialDff')
    flData = trialDffTuningShifted;
elseif strcmp(flType, 'expDff')
    flData = expDffTuningShifted;
end
tuningFl = smoothdata(flData, 1, 'gaussian', smWin, 'omitnan');

% Find max and min for each trial and EB wedge
minVals = [];
maxVals = [];
for iTrial = 1:nTrials
    minVals(:, iTrial) = min(tuningFl(:, :, iTrial), [], 1, 'omitnan'); % --> [wedge, trial]
    maxVals(:, iTrial) = max(tuningFl(:, :, iTrial), [], 1, 'omitnan'); % --> [wedge, trial]
end
tuningAmp = maxVals - minVals; % --> [wedge, trial]

% Group into pre- and post- drug trials
titleStr = [];
if ~isempty(drugTrials)
    blockAmps = [];
    for iTrial = 1:numel(drugTrials)
        preStim = drugTrials(iTrial) - 1;
        postStim = drugTrials(iTrial);
        if ~currExpData.usingPanels(postStim)
            postStim = postStim + 1;
            titleStr{iTrial} = 'In darkness';
        else
            titleStr{iTrial} = 'With visual stim';
        end
        blockAmps(:, :, iTrial) = tuningAmp(:, [preStim, postStim]); % --> [wedge, trial, stim trial num]
    end
    blockAmps = permute(blockAmps, [2, 1 3]); % --> [trial, wedge, stim trial Num]
end

% Create figure;
f = figure(6); clf;
f.Color = [1 1 1];
if ~isempty(figPos)
    f.Position(3:4) = figPos;
    if sum(f.Position([2, 4])) > 1080
        f.Position(2) = 50;
    end
end

nPlots = size(flData, 2);
for iPlot = 1:nPlots
    ax = subaxis(nPlots, 1, iPlot, 'mt', MT, 'mb', MB, 'sv', SV, 'mr', MR, ...
        'ml', ML, 'sh', SH);
    hold on;
    
    if strcmp(plotType, 'amplitude')
        plot(1:size(tuningAmp, 2), tuningAmp(iPlot, :), 'o', 'color', 'k', 'markerfacecolor', 'k')
        ylim([0, 1.2 * max(tuningAmp(:))]);
    elseif strcmp(plotType, 'minmax')
        for iTrial = 1:nTrials
            plot([1, 1] .* iTrial, [minVals(iPlot, iTrial), maxVals(iPlot, iTrial)], 'k')
        end
        clear plt;
        plt(1) = plot(1:size(minVals, 2), minVals(iPlot, :), 'o', 'color', 'b', ...
                'markerfacecolor', 'b');
        plt(2) = plot(1:size(maxVals, 2), maxVals(iPlot, :), 'o', 'color', 'm', ...
                'markerfacecolor', 'm');
        ylim([min(minVals(:)) - abs(min(minVals(:))), max(maxVals(:))] .* [1, 1.2]) 
        if iPlot == 1
            legend(plt, {'Min', 'Max'}, 'autoupdate', 'off', 'fontsize', 10)
        end
    end
    
    % Indicate transitions to drug application trials
    if ~isempty(drugTrials)
        for iTrial = 1:numel(drugTrials)
            yL = ylim();
            if strcmp(currExpData.visualStim{drugTrials(iTrial)}, 'yes')
                lineColor = 'g';
            else
                lineColor = rgb('yellow');
            end
            plot([-0.5, -0.5] + drugTrials(iTrial), yL, 'color', lineColor, 'linewidth', 3)
        end
    end
    
    % Adjust axes
    if iPlot == nPlots
        xlabel('Trial number', 'fontsize', 13)
    else
        ax.XTick = [];
        ax.XTickLabel = [];
    end
    ax.FontSize = 12;
    xL = xlim;
    xlim([0.5, size(tuningAmp, 2) + 0.5])
%     box off;
end

figTitleStr = {expID, ['EB wedge ', flType, ' tuning curve amplitudes']};
if ~isempty(drugTrials)
    figTitleStr = [figTitleStr, ['(Green lines = ', currExpData.drugName{drugTrials(1)}, ...
            ' application)']];
end
h = suptitle(figTitleStr);
h.FontSize = 16;


% Also create an amplitude summary figure
nStimTrials = size(blockAmps, 3);
g = figure(7);clf;
g.Color = [1 1 1];
g.Position = [g.Position(1:2), nStimTrials * 300, 450];
clear allAx;
for iTrial = 1:nStimTrials
    subaxis(1, nStimTrials, iTrial, 'mr', 0.04, 'ml', 72 / f.Position(3), 'mt', 0.15, 'mb', 0.11, ...
            'sh', 0.08); 
    hold on; box off
    plot(blockAmps(:, :, iTrial), 'o-')
    ax = gca;
    ax.XTick = [1 2];
    ax.XTickLabel = {'Before', 'After'};
    ax.FontSize = 14;
    xlim([0.75 2.25]);
    title(titleStr{iTrial}, 'fontsize', 12, 'verticalalignment', 'top');
    
    % Plot mean amplitude for each group
    plot([mean(blockAmps(1, :, iTrial)), mean(blockAmps(2, :, iTrial))], ...
            'o-', 'Color', 'k', 'linewidth', 5)
    
    % Track Y limits of each axis
    if iTrial == 1
       figYLims = ylim();
       ylabel('Max - min tuning curve dF/F', 'fontsize', 14) 
    else
       currYLims = ylim();
       figYLims = [min([figYLims(1), currYLims(1)]), max([figYLims(2), currYLims(2)])];
    end
    allAx(iTrial) = ax;
end

% Match Ylims across figures
for iAx = 1:nStimTrials
    allAx(iAx).YLim = figYLims;
end

% Add title
suptitle({expID, ['Tuning curve amplitudes before and after ', ...
        currExpData.drugName{drugTrials(1)}]})

% Save figure
if saveFig
    figTitle = [expID, '_visual_tuning_minmax'];
    save_figure(f, figDir, figTitle);
    figTitle = [expID, '_visual_tuning_amplitude_summary'];
    save_figure(g, figDir, figTitle);
end

catch ME; rethrow(ME); end

%% PLOT ALL TUNING CURVES ON OVERLAID POLAR PLOTS

saveFig = 0;

smWin = 5;

flType = 'expDff';

matchRLims = 1;
figPos = [];
% figPos = [1850 320];

% Adjust plot spacing and margins
SV = 0.03;
SH = 0.05;
ML = 0.03;
MR = 0.03;
MT = 0.3;
MB = 0.05;

try 
    
% Get correct source data
if strcmp(flType, 'rawFl')
    flData = rawFlTuningShifted;
elseif strcmp(flType, 'trialDff')
    flData = trialDffTuningShifted;
elseif strcmp(flType, 'expDff')
    flData = expDffTuningShifted;
end
tuningFl = smoothdata(flData, 1, 'gaussian', smWin, 'omitnan');

% Create figure;
f = figure(8); clf;
f.Color = [1 1 1];
if ~isempty(figPos)
    f.Position(3:4) = figPos;
    if sum(f.Position([2, 4])) > 1080
        f.Position(2) = 50;
    end
end

nPlots = nTrials;
for iPlot = 1:nPlots
    
    % Create axes
    ax = subaxis(1, nPlots, iPlot, 'mt', MT, 'mb', MB, 'sv', SV, 'mr', MR, ...
        'ml', ML, 'sh', SH);
    hold on;
    plotX = deg2rad(-180:3.75:(180 - 3.75));
    pax = polaraxes(); hold on
    pax.Position = ax.Position;
    ax.Visible = 'off';
    
    % Plot data
    currData = tuningFl(:, :, iPlot);
    cm = hsv(nRois) .* 0.9;
    for iWedge = 1:nRois
%         plotData = currData(:, iWedge) - min(currData(:));
        plotData = currData(:, iWedge);
        polarplot(plotX, plotData, 'color', cm(iWedge, :), 'linewidth', 1.5);
    end
    
    % Format axes
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.ThetaTick = [0:45:179, 225:45:359];
    pax.ThetaTickLabel = {'0' '+45', '+90' '+135' '-135', '-90', '-45'};
    pax.RTick = [];
    pax.FontSize = 12;
    if matchRLims
        pax.RLim(2) = 1.01 * max(as_vector(smoothdata(flData, 1, 'gaussian', smWin)));
    end
    title(['#', num2str(iPlot)], 'fontSize', 14)
    
    % Indicate transitions to drug application trials
    if ~isempty(drugTrials) && ismember(iPlot, drugTrials)
        pax.Title.String = [pax.Title.String, ' (post-', currExpData.drugName{iPlot}, ')'];
    end
end

titleStr = [expID, '  —  EB wedge ', flType, ' tuning curves'];
if matchRLims
    titleStr = [titleStr, '  (axis limits matched across plots)'];
else
    titleStr = [titleStr, '  (axis limits set independently across plots)'];
end
h = suptitle(titleStr);
h.FontSize = 16;

% Save figure
if saveFig
    figTitle = [expID, '_visual_tuning_polarPlots'];
    save_figure(f, figDir, figTitle);
end

catch ME; rethrow(ME); end



%% PLOT WHOLE-TRIAL VECTOR STRENGTH FOR EACH EB WEDGE

smWin = 5;

flType = 'expDff';

figPos = [550 950];

% Adjust plot spacing and margins
SV = 0.03;
SH = 0.03;
ML = 0.08;
MR = 0.02;
MT = 0.05;
MB = 0.07;

% Create figure;
f = figure(9); clf;
f.Color = [1 1 1];
if ~isempty(figPos)
    f.Position(3:4) = figPos;
    if sum(f.Position([2, 4])) > 1080
        f.Position(2) = 50;
    end
end

nPlots = size(fullTrialVectorStrength, 2);
for iPlot = 1:nPlots
    ax = subaxis(nPlots, 1, iPlot, 'mt', MT, 'mb', MB, 'sv', SV, 'mr', MR, ...
        'ml', ML, 'sh', SH);
    hold on;
    
    plot(1:size(fullTrialVectorStrength, 1), fullTrialVectorStrength(:, iPlot), 'o', 'color', 'k', ...
            'markerfacecolor', 'k')
    ylim([0, 1.2 * max(fullTrialVectorStrength(:))]);
    
    % Indicate transitions to drug application trials
    if ~isempty(drugTrials)
        for iTrial = 1:numel(drugTrials)
            yL = ylim();
           plot([-0.5, -0.5] + drugTrials(iTrial), yL, 'color', 'g', 'linewidth', 3)
        end
    end
    
    % Adjust axes
    if iPlot == nPlots
        xlabel('Trial number', 'fontsize', 13)
    else
        ax.XTick = [];
        ax.XTickLabel = [];
    end
    ax.FontSize = 12;
    xL = xlim;
    xlim([0.5, size(fullTrialVectorStrength, 1) + 0.5])
%     box off;
end


% f = figure(8); clf;
% f.Color = [1 1 1];
% 
% % Adjust plot spacing and margins
% SV = 0.06;
% SH = 0.04;
% ML = 0.04;
% MR = 0.04;
% MT = 0.05;
% MB = 0.07;
% count = 1;
% for iTrial = 1:nTrials
%     for iRoi = 1:nRois
%         clear pax;
%         ax = subaxis(nTrials, nRois, count, 'mt', MT, 'mb', MB, 'sv', SV, 'mr', MR, ...
%                 'ml', ML, 'sh', SH);
%         count = count + 1;
%         pax = polaraxes();
%         pax.Position = ax.Position;
% %         pax.Position = pax.Position - [0.04 0.025 0 0];
% %         pax.Position = pax.Position .* [1 1 1.2 1.2];
%         ax.Visible = 'off';
%         
%         polarplot([1, 1] .* (fullTrialVectorPhase(iTrial, iRoi)), [0, ...
%                 fullTrialVectorStrength(iTrial, iRoi)], 'color', 'k', 'linewidth', 4)
%         pax.ThetaZeroLocation = 'top';
% % %         pax.ThetaDir = 'clockwise';
%         pax.ThetaTick = [0:45:179, 225:45:359];
% %         pax.ThetaTickLabel = {'0' '+45', '+90' '+135' '-135', '-90', '-45'};
%         pax.RTick = [];
%     end
% end






























