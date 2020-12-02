startDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';
parentDir = uigetdir(startDir, 'Select an experiment directory');

analysisDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\Analysis';


load(fullfile(parentDir, 'analysis_data.mat'));


% Add optic flow info to mD struct
if exist(fullfile(parentDir, 'flowMags.mat'), 'file')
    load(fullfile(parentDir, 'flowMags.mat'));
    for iTrial = 1:numel(mD)
        mD(iTrial).flowData = meanFlowMags{iTrial};
        mD(iTrial).flowFrameDur = median(diff(mD(iTrial).ftFrameTimes));
        mD(iTrial).flowFrameTimes = (1:1:numel(meanFlowMags{iTrial})) * mD(iTrial).flowFrameDur;
    end
end

%% COMBINE DATA FROM A BLOCK OF COMPATIBLE TRIALS

blTrials = [1:15];


flowSmWin = 30;
moveThresh = 0.05;
moveDistThresh = 2;


sectionBreakTrials = [6 10 13 14];
sectionColorLabels = {'500 uM DA', 'Stopped chilling saline', 'No visual stim', 'Back to visual stim'};%{'100 uM DA', '250 uM DA', '500 uM DA', 'No visual stim', 'Back to visual stim'};
sectionBreakColors = [1:numel(sectionBreakTrials)];
fullColorList = [rgb('magenta'); rgb('green'); rgb('gold'); rgb('cyan') * 0.8; rgb('red') * 0.8];
sectionBreakColorList = fullColorList(1:numel(sectionBreakTrials), :);


expType = 'PB';%'EB';%

bl = extract_block_data(mD, blTrials, ...
        'flowSmWin', flowSmWin, 'moveThresh', moveThresh, 'ExpType', expType);

    
%% PLOT SUMMARY OF MOVEMENT THROUGHOUT EXPERIMENT

saveFig = 0;

try
    
% Heatmap of slightly smoothed flow throughout experiment
f = figure(1);clf;
f.Color = [1 1 1];
ax = gca;
allFlow = bl.meanVolFlow;
allFlow([1:5, size(allFlow, 1)-5:size(allFlow, 1)], :) = nan;
allFlowSm = smoothdata(allFlow, 1, 'gaussian', 4, 'omitnan');
allFlowSm = repeat_smooth(allFlow, 20, 'dim', 1, 'smwin', 5);
imagesc([0, bl.trialDuration], [1, numel(bl.trialNum)], allFlowSm')
title('Optic flow')
ax.YTick = 1:numel(bl.trialNum);
xlabel('Time (sec)')
ylabel('Trial')
ax.FontSize = 12;

% Plot avg flow value and percentage of movement vols for each trial
f = figure(2); clf; hold on;
f.Color = [1 1 1];
f.Position(3:4) = [1200 400];
ax = subaxis(1, 1, 1, 'mt', 0.1, 'mb', 0.16, 'ml', 0.08, 'mr', 0.08); hold on
avgTrialFlow = mean(allFlow, 1, 'omitnan');
plotX = 1:numel(avgTrialFlow);
plot(plotX, avgTrialFlow, '-o', 'color', 'b', 'linewidth', 1, 'markersize', 8)
ax.YColor = [0 0 1];
stimTrials = find([bl.usingOptoStim]);
for iStim = 1:numel(stimTrials)
    currTrialInd = stimTrials(iStim);
    if bl.usingPanels(currTrialInd)
        shadeColor = [0 1 0];
    else
        shadeColor = [1 0 0];
    end
    plot_stim_shading([currTrialInd - 0.5, currTrialInd + 0.5], 'color', shadeColor);
end

% Indicate stim trials
optoStimTrials = find([bl.usingOptoStim]);
if ~isempty(optoStimTrials)
    sectionBreakColorList = [rgb('green'); rgb('red')];
    sectionColorLabels = {'Visual+opto', 'Opto-only'};
    sectionBreakColors = [];
    sectionBreakTrials = [];
    for iTrial = 1:numel(optoStimTrials)
        sectionBreakTrials(end + 1) = optoStimTrials(iTrial);
        if bl.usingPanels(optoStimTrials(iTrial))
            sectionBreakColors(end + 1) = 1;
        else
            sectionBreakColors(end + 1) = 2;
        end
    end%iTrial
end%if
    
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

ax.FontSize = 12;
xlabel('Trial', 'fontsize', 16)
xlim([0 plotX(end) + 1])
ax.XTick = 1:numel(bl.trialNum);
ax.XTickLabel = [bl.trialNum];
smFlow = repeat_smooth(allFlow, 20, 'dim', 1, 'smwin', flowSmWin);
smFlow = smFlow - min(smFlow(:));
moveVols = smFlow > moveThresh;
trialMovePercent = (sum(moveVols) ./ size(moveVols, 1)) * 100;
ylabel('Mean optic flow (AU)', 'fontsize', 16)

% Plot proportion of volumes when fly was moving in each trial
yyaxis('right');
ax.YColor = [1 0 0];
ylabel('% movement', 'fontsize', 16)
plot(trialMovePercent, '-*', 'color', 'r', 'linewidth', 1, 'markersize', 10)
ylim([0 100])
yyaxis('left');
ax.YColor = [0 0 1];

% Add title
figTitleStr = [bl.expID, '  -  Fly movement summary (green = opto + visual, red = opto only)'];
t = title(figTitleStr);
t.Position(2) = t.Position(2) * 1.02;

% Save figure
saveFileName = 'fly_movement_summary';
if saveFig
   saveDir = fullfile(analysisDir, bl.expID);
   save_figure(f, saveDir, saveFileName);
end

catch ME; rethrow(ME); end

%% PLOT SINGLE-TRIAL HEATMAPS FOR ENTIRE BLOCK

saveFig =  0;

omitMoveVols = 0;

% sourceData = bl.wedgeRawFlArr;
sourceData = bl.wedgeDffArr;
% % sourceData = bl.wedgeZscoreArr;
sourceData = bl.wedgeExpDffArr;

% sourceData =  bl.dffArr;
% sourceData =  bl.rawFlArr;
% sourceData =  bl.zscoreArr;
% sourceData = bl.expDffArr;

figPos = [650 980];

try 
    
% Generate figure labels and save file name
if isequal(sourceData, bl.wedgeDffArr) || isequal(sourceData, bl.dffArr)
    figTitleText = 'dF/F';
    saveFileName = 'single_trial_heatmaps_dff';
elseif isequal(sourceData, bl.wedgeRawFlArr) || isequal(sourceData, bl.rawFlArr)
    figTitleText = 'raw F';
    saveFileName = 'single_trial_heatmaps_rawF';
elseif isequal(sourceData, bl.wedgeZscoreArr) || isequal(sourceData, bl.zscoreArr)
    figTitleText = 'Z-scored dF/F';
    saveFileName = 'single_trial_heatmaps_zscored-dff';
elseif isequal(sourceData, bl.wedgeExpDffArr)  || isequal(sourceData, bl.expDffArr)
    figTitleText = 'full exp dF/F';
    saveFileName = 'single_trial_heatmaps_exp-dff';
else
%     warndlg('Warning: sourceData mismatch');
    figTitleText = '';
    saveFileName = 'single_trial_heatmaps_unknown';
end

% Replace volumes when fly was moving with nan if necessary
if omitMoveVols
   moveDistThreshVols = round(moveDistThresh * bl.volumeRate);
   moveDistArr = permute(repmat(bl.moveDistVols, 1, 1, size(sourceData, 2)), [1 3 2]);
   sourceData(moveDistArr < moveDistThreshVols) = nan; 
   saveFileName = [saveFileName, '_no-movement'];
end

nTrials = numel([bl.trialNum]);

% To give all figures the same color scale
sourceData(1, 1, :) = min(as_vector(sourceData(:, :, ~[bl.usingOptoStim])), [], 'omitnan');
sourceData(end, end, :) = max(sourceData(:), [], 'omitnan');

f = figure(10);clf;
f.Color = [1 1 1];
if ~isempty(figPos)
    f.Position(3:4) = figPos;
end
colormap('parula')
for iTrial = 1:nTrials 
    
   subaxis(nTrials, 4, [1:3] + (4 * (iTrial - 1)), 'mt', 0, 'mb', 0.06, 'sv', 0.001, 'mr', 0.07);
   
   imagesc([0, bl.trialDuration], [1, 8], smoothdata(sourceData(:, :, iTrial), 1, 'gaussian', 3)');
   hold on
   invertedPhaseData = (-bl.dffVectAvgRad/pi) * 4 + 4.5;
   plot(bl.volTimes, invertedPhaseData(:, iTrial), 'color', 'r')
   
   % Fill in opto stim trials with solid color
   ax = gca;
   ax.FontSize = 12;
   if bl.usingOptoStim(iTrial) && bl.usingPanels(iTrial)
       colormap(ax, [0 1 0])
   elseif bl.usingOptoStim(iTrial)
       colormap(ax, [1 0 0])
   end
   
   % Label X axis on final trial
   if iTrial < nTrials
       ax.XTickLabel = [];
   else
       xlabel('Time (sec)', 'FontSize', 12);
   end
   
   % Label each plot with trial number
   ax.YTickLabel = [];
   t = ylabel(num2str(bl.trialNum(iTrial)), 'rotation', 0, 'FontSize', 12);
   t.HorizontalAlignment = 'right';
   t.VerticalAlignment = 'middle';
   t.Position(1) = t.Position(1) * 2;
   
   
   % Add fly movement summary plot to the side of each trial
   ax = subaxis(nTrials, 4, 4*iTrial);
   normTrialFlow = avgTrialFlow ./ max(avgTrialFlow);
%    b = barh([normTrialFlow(iTrial), trialMovePercent(iTrial) / 100; nan, nan]);
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
   end
   
end

% Add Y-axis label to entire figure
tAx = axes(f, 'Position', [0 0 1 1]);
t = text(tAx, 0.03, 0.5, 'Trial', 'units', 'normalized');
t.Rotation = 90;
t.FontSize = 14;
tAx.Visible = 'off';

% Add figure title
figTitleStr = {[bl.expID, '  -  single trial ', figTitleText, ' across EB wedges'], ...
        'green = opto + visual, red = opto only'};
if omitMoveVols
    figTitleStr = [figTitleStr(1), {['Excl. volumes within ', num2str(moveDistThresh), ...
            ' sec of movement']}, figTitleStr{2}];
end
h = suptitle(figTitleStr);
h.FontSize = 14;

% Save figure
if saveFig
   saveDir = fullfile(analysisDir, bl.expID);
   save_figure(f, saveDir, saveFileName);
end

catch ME; rethrow(ME); end

%% ===================================================================================================
%% Plot tuning heatmaps for each wedge across trials
%===================================================================================================
% opts = [];

saveFig = 1;

omitMoveVols = 0;

showStimTrials = 0;

sourceData = bl.wedgeRawFlArr;
% sourceData = bl.wedgeDffArr;
% sourceData = bl.wedgeZscoreArr;
sourceData = bl.wedgeExpDffArr;

try
    
% Generate figure labels and save file name
if isequal(sourceData, bl.wedgeDffArr)
    figTitleText = 'dF/F';
    saveFileName = '2D_tuning_dff';
elseif isequal(sourceData, bl.wedgeRawFlArr)
    figTitleText = 'raw F';
    saveFileName = '2D_tuning_rawF';
elseif isequal(sourceData, bl.wedgeZscoreArr)
    figTitleText = 'Z-scored dF/F';
    saveFileName = '2D_tuning_zscored-dff';
elseif isequal(sourceData, bl.wedgeExpDffArr)
    figTitleText = 'full exp dF/F';
    saveFileName = '2D_tuning_exp-dff';
else
    warndlg('Warning: sourceData mismatch');
    figTitleText = '';
    saveFileName = '2D_tuning_unknown';
end

% Replace volumes when fly was moving with nan if necessary
if omitMoveVols
   moveDistThreshVols = round(moveDistThresh * bl.volumeRate);
   moveDistArr = permute(repmat(bl.moveDistVols, 1, 1, 8), [1 3 2]);
   sourceData(moveDistArr < moveDistThreshVols) = nan; 
   saveFileName = [saveFileName, '_no-movement'];
end

nTrials = numel(bl);

% Get mean panels pos data
panelsPosVols = [];
for iVol = 1:size(sourceData, 1)
    [~, currVol] = min(abs(bl.panelsFrameTimes - bl.volTimes(iVol)));
    panelsPosVols(iVol) = bl.panelsPosX(currVol);
end
meanData = [];
for iPos = 1:numel(unique(bl.panelsPosX))
    meanData(iPos, :, :) = ...
            mean(sourceData(panelsPosVols == (iPos - 1), :, :), 1, 'omitnan'); % --> [barPos, wedge, trial]    
end
meanDataSmooth = smoothdata(meanData, 1, 'gaussian', 3, 'omitnan');

% Determine subplot grid size
nPlots = size(meanDataSmooth, 2);
if nPlots == 3
    plotPos = [2 2]; % Because [1 3] looks bad
else
    plotPos = numSubplots(nPlots);
end

% Create figure and plots
f = figure(1);clf;
f.Color = [1 1 1];
for iPlot = 1:nPlots
    subaxis(plotPos(1), plotPos(2), iPlot, 'ml', 0.05, 'mr', 0.05);
    
    currData = squeeze(meanDataSmooth(:, iPlot, :)); % --> [barPos, trial]
    if ~showStimTrials
        currData = currData(:, ~[bl.usingOptoStim]);
    end
    currDataShift = [currData(92:96, :); currData(1:91, :)];
    
    % Plot data
    plotX = -180:3.75:(180 - 3.75);
    imagesc(plotX, [1 size(currData, 2)], ...
            smoothdata(currDataShift, 1, 'gaussian', 3)');
    hold on; plot([0 0], [0 size(currData, 2) + 1], 'color', 'k', 'linewidth', 2)
    ylim([0.5 size(currData, 2) + 0.5])
%     colorbar
    xlabel('Bar position (degrees from front of fly)');
    ylabel('Trial');
    
    if showStimTrials
        % Label opto stim trials
        optoStimTrials = find([bl.usingOptoStim]);
        for iTrial = 1:numel(optoStimTrials)
            if bl.usingPanels(optoStimTrials(iTrial))
                lineColor = 'green';
            else
                lineColor = 'red';
            end
            xL = xlim() + [2 -2];
            plotX = [xL(1), xL(2), xL(2), xL(1), xL(1)];
            plotY = [-0.5 -0.5 0.5 0.5 -0.5] + optoStimTrials(iTrial);
            plot(plotX, plotY, 'color', lineColor, 'linewidth', 1.5)
            
        end
        
        % Add Y tick labels
        ax = gca;
        ax.YTick = 1:size(currData, 2);
        ax.YTickLabels = blTrials;
        
    else
        % Indicate missing opto stim trials
            skippedTrials = 0;
            optoStimTrials = find([bl.usingOptoStim]);
            for iTrial = 1:numel(optoStimTrials)
                plotY = optoStimTrials(iTrial) - skippedTrials - 0.5;
                skippedTrials = skippedTrials + 1;
                if bl.usingPanels(optoStimTrials(iTrial))
                    lineColor = 'green';
                else
                    lineColor = 'red';
                end
                xL = xlim;
                plot(xL, [plotY, plotY], 'color', lineColor, 'linewidth', 2);
            end
            
            % Skip omitted trial numbers in Y tick labels
            ax = gca();
            ax.YTick = 1:numel(bl.usingOptoStim);
            yTickLabel = blTrials;
            ax.YTickLabel = yTickLabel(~bl.usingOptoStim);
            
    end%if
end%iPlot

% Add figure title
figTitleStr = {[bl.expID, '  -  Visual tuning of EB wedges (mean ', figTitleText, ')'], ...
        'Green line  =  trial with opto + visual stim', ...
        '  Red line  = trial with opto stim only'};
if omitMoveVols
    figTitleStr{1} = [figTitleStr{1}, '  -  excluding volumes within ', num2str(moveDistThresh), ...
            ' sec of movement'];
end
h = suptitle(figTitleStr);
h.FontSize = 14;
    
% Save figure
if saveFig
   saveDir = fullfile(analysisDir, bl.expID);
   save_figure(f, saveDir, saveFileName);
end

catch ME; rethrow(ME); end

%% ===================================================================================================    
%% Plot as lines instead of using imagesc
% ===================================================================================================

saveFig = 1;

omitMoveVols = 0;

smWin = 3;

sourceData = bl.wedgeRawFlArr;
sourceData = bl.wedgeDffArr;
% sourceData = bl.wedgeZscoreArr;
sourceData = bl.wedgeExpDffArr;


figSize = [];
figSize = [1250 980];

try
    
% Generate figure labels and save file name
if isequal(sourceData, bl.wedgeDffArr)
    figTitleText = 'dF/F';
    saveFileName = 'tuning_curves_dff';
elseif isequal(sourceData, bl.wedgeRawFlArr)
    figTitleText = 'raw F';
    saveFileName = 'tuning_curves_rawF';
elseif isequal(sourceData, bl.wedgeZscoreArr)
    figTitleText = 'Z-scored dF/F';
    saveFileName = 'tuning_curves_zscored-dff';
elseif isequal(sourceData, bl.wedgeExpDffArr)
    figTitleText = 'full exp dF/F';
    saveFileName = 'tuning_curves_exp-dff';
else
    warndlg('Warning: sourceData mismatch');
    figTitleText = '';
    saveFileName = 'tuning_curves_unknown';
end

% Replace volumes when fly was moving with nan if necessary
if omitMoveVols
   moveDistThreshVols = round(moveDistThresh * bl.volumeRate);
   moveDistArr = permute(repmat(bl.moveDistVols, 1, 1, 8), [1 3 2]);
   sourceData(moveDistArr < moveDistThreshVols) = nan; 
   saveFileName = [saveFileName, '_no-movement'];
end

% Get mean panels pos data
panelsPosVols = [];
for iVol = 1:size(sourceData, 1)
    [~, currVol] = min(abs(bl.panelsFrameTimes - bl.volTimes(iVol)));
    panelsPosVols(iVol) = bl.panelsPosX(currVol);
end
plotData = [];
for iPos = 1:numel(unique(bl.panelsPosX))
    plotData(iPos, :, :) = ...
            mean(sourceData(panelsPosVols == (iPos - 1), :, :), 1, 'omitnan'); % --> [barPos, wedge, trial]    
end

% Shift data so center of plot is directly in front of the fly
shiftData = cat(1, plotData(92:96, :, :), plotData(1:91, :, :));
shiftDataSm = repeat_smooth(shiftData, 10, 'smWin', smWin);

% Offset data so that each plot is centered at zero
shiftDataOffset = [];
for iTrial = 1:size(shiftDataSm, 3)
   for iWedge = 1:size(shiftDataSm, 2)
      currData = shiftDataSm(:, iWedge, iTrial); % --> [barPos]
      shiftDataOffset(:, iWedge, iTrial) = currData - mean(currData, 'omitnan');
   end    
end

% Find max and min values 
% yMax = max(shiftDataOffset(800:end));
% yMin = min(shiftDataOffset(800:end));
yMax = max(shiftDataOffset(:), [], 'omitnan');
yMin = min(shiftDataOffset(:), [], 'omitnan');
range = max(abs([yMax, yMin]), [], 'omitnan') *  2.5;

% Create figure and plots
f = figure(11);clf;
f.Color = [1 1 1];
if ~isempty(figSize)
   f.Position(3:4) = figSize; 
end
allAx = {}; plotLines = {};
for iPlot = 1:8
    allAx{iPlot} = subaxis(1, nPlots, iPlot, 'mt', 0, 'ml', 0.05, 'mr', 0.03, 'mb', 0.08, 'sh', 0.05);
    hold on;
    currData = squeeze(shiftDataOffset(:, iPlot, :)); % --> [barPos, trial]
    currData = currData(:, ~[bl.usingOptoStim]);
    
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
    allAx{iPlot}.YTickLabel = yTickLabel(~fliplr(bl.usingOptoStim));
    
    % Indicate missing opto stim trials
    skippedTrials = 0;
    optoStimTrials = find([bl.usingOptoStim]);
    if ~isempty(optoStimTrials)
        sectionBreakColorList = [rgb('green'); rgb('red')];
        sectionColorLabels = {'Visual+opto', 'Opto-only'};
        sectionBreakColors = [];
        sectionBreakTrials = [];
        for iTrial = 1:numel(optoStimTrials)
            sectionBreakTrials(end + 1) = optoStimTrials(iTrial) - skippedTrials;
            skippedTrials = skippedTrials + 1;
            if bl.usingPanels(optoStimTrials(iTrial))
                sectionBreakColors(end + 1) = 1;
            else
                sectionBreakColors(end + 1) = 2;
            end
        end%iTrial
    end%if
    
    % Plot lines dividing sections if necessary
    yL = ylim();
    if ~isempty(sectionBreakTrials)
        offsetsRev = offsets(end:-1:1);
        for iColor = 1:size(sectionBreakColorList, 1)
            currColor = sectionBreakColorList(iColor, :);
            currSectionBreaks = sectionBreakTrials(sectionBreakColors == iColor);
            for iSection = 1:numel(currSectionBreaks)
                currTrialNum = currSectionBreaks(iSection);
                plotOffset = offsetsRev(currTrialNum - 1) - (range/2);
                xx = [plotX(1), plotX(end)];
                yy = [plotOffset plotOffset];
                if iColor == 1 && iSection == 1
                    plotLines{iPlot} = plot(xx, yy, 'color', currColor, 'linewidth', 2);
                elseif iSection == 1
                    plotLines{iPlot}(end + 1) = plot(xx, yy, 'color', currColor, 'linewidth', 2);
                else
                    plot(xx, yy, 'color', currColor, 'linewidth', 2);
                end
            end%iSection
        end%iColor 
    end%if  
end%iPlot

% Add figure title
figTitleStr = [bl.expID, '  -  Visual tuning of EB wedges (mean ', figTitleText, ')'];    
if omitMoveVols
    figTitleStr{1} = [figTitleStr{1}, '  -  excluding volumes within ', num2str(moveDistThresh), ...
            ' sec of movement'];
end
h = suptitle(figTitleStr);
h.FontSize = 14;

% Add legend for any opto stim trials or section breaks
for iPlot = 1:8
    clear lgd;
    lgd = legend(allAx{iPlot}, plotLines{iPlot}, sectionColorLabels, 'location', 'none');
    lgd.FontSize = 12;
    lgd.Orientation = 'horizontal';
    lgd.Box = 'off'; 
    drawnow()
    lgd.Position(3) = lgd.Position(3) * 1.2;
    axPos = allAx{1}.Position;
    lgd.Position(1) = axPos(1);
    lgd.Position(2) = sum(axPos([2 4])) * 1.005;
    if iPlot < 8
        lgd.Visible = 'off';
    end
end

% Save figure
if saveFig
   saveDir = fullfile(analysisDir, bl.expID);
   save_figure(f, saveDir, saveFileName);
end

catch ME; rethrow(ME); end

%% ===================================================================================================
%% Plot min and max values from the tuning curves
%===================================================================================================

saveFig = 0;

omitMoveVols = 0;

smWin = 3;

sourceData = bl.wedgeRawFlArr;
sourceData = bl.wedgeDffArr;
% sourceData = bl.wedgeZscoreArr;
sourceData = bl.wedgeExpDffArr;

figSize = [];
figSize = [400 925];

try
    
% Generate figure labels and save file name
if isequal(sourceData, bl.wedgeDffArr)
    figTitleText = 'dF/F';
    saveFileName = 'tuning_curve_amplitudes_dff';
elseif isequal(sourceData, bl.wedgeRawFlArr)
    figTitleText = 'Raw F';
    saveFileName = 'tuning_curve_amplitudes_rawF';
elseif isequal(sourceData, bl.wedgeZscoreArr)
    figTitleText = 'Z-scored dF/F';
    saveFileName = 'tuning_curve_amplitudes_zscored-dff';
elseif isequal(sourceData, bl.wedgeExpDffArr)
    figTitleText = 'full exp dF/F';
    saveFileName = 'tuning_curve_amplitudes_exp-dff';
else
    warndlg('Warning: sourceData mismatch');
    figTitleText = '';
    saveFileName = 'tuning_curve_amplitudes_unknown';
end

% Replace volumes when fly was moving with nan if necessary
if omitMoveVols
   moveDistThreshVols = round(moveDistThresh * bl.volumeRate);
   moveDistArr = permute(repmat(bl.moveDistVols, 1, 1, 8), [1 3 2]);
   sourceData(moveDistArr < moveDistThreshVols) = nan; 
   saveFileName = [saveFileName, '_no-movement'];
end

% Get mean panels pos data
panelsPosVols = [];
for iVol = 1:size(bl.wedgeDffArr, 1)
    [~, currVol] = min(abs(bl.panelsFrameTimes - bl.volTimes(iVol)), [], 'omitnan');
    panelsPosVols(iVol) = bl.panelsPosX(currVol);
end
tuningData = [];
for iPos = 1:numel(unique(bl.panelsPosX))
    tuningData(iPos, :, :) = ...
            mean(sourceData(panelsPosVols == (iPos - 1), :, :), 1, 'omitnan'); % --> [barPos, wedge, trial]    
end

% Smooth data, then find the max and min for each trial and wedge
tuningDataSm = repeat_smooth(tuningData, 10, 'smWin', smWin);
minVals = []; maxVals = [];
for iTrial = 1:size(tuningDataSm, 3)
    minVals(:, iTrial) = min(tuningDataSm(:, :, iTrial), [], 1, 'omitnan'); % --> [wedge, trial]
    maxVals(:, iTrial) = max(tuningDataSm(:, :, iTrial), [], 1, 'omitnan'); % --> [wedge, trial]
end
xTickLabel = 1:size(minVals, 2);

% Cut out opto stim trials
optoStimTrials = [bl.usingOptoStim];
postStimTrials = [0, optoStimTrials(1:end - 1)];
postStimTrials = postStimTrials(~optoStimTrials);
optoStimPanels = bl.usingPanels;
prevStimPanels = [0, optoStimPanels(1:end - 1)];
prevStimPanels = prevStimPanels(~optoStimTrials);
xTickLabel = xTickLabel(~optoStimTrials);
minVals = minVals(:, ~optoStimTrials);  % --> [wedge, trial]
maxVals = maxVals(:, ~optoStimTrials);  % --> [wedge, trial]
tuningAmp = maxVals - minVals;          % --> [wedge, trial]

% Plot change in min and max over time, omitting opto stim trials
f = figure(2);clf;
f.Color = [1 1 1];
if ~isempty(figSize)
    f.Position(2) = 50;
   f.Position(3:4) = figSize; 
end
globalMax = 0;
nPlots = size(sourceData, 2);
for iPlot = 1:nPlots
    ax = subaxis(nPlots, 1, iPlot, 'mt', 0.08, 'mb', 0.06, 'sv', 0.03, 'mr', 0.03, 'ml', 0.08);
    hold on;
%     plot(1:size(minVals, 2), minVals(iPlot, :), 'o', 'color', 'b')
%     plot(1:size(maxVals, 2), maxVals(iPlot, :), 'o', 'color', 'm')
    plot(1:size(tuningAmp, 2), tuningAmp(iPlot, :), 'o', 'color', 'k')
    
    % Indicate stim trials
    stimTrials = find(postStimTrials);
    skipTrials = 0;
    if ~isempty(stimTrials)
        sectionBreakColorList = [rgb('green'); rgb('red')];
        sectionColorLabels = {'Visual+opto', 'Opto-only'};
        sectionBreakColors = [];
        sectionBreakTrials = [];
        for iTrial = 1:numel(stimTrials)
            sectionBreakTrials(end + 1) = stimTrials(iTrial) - skipTrials;
            skipTrials = skipTrials + 1;
            if bl.usingPanels(stimTrials(iTrial))
                sectionBreakColors(end + 1) = 1;
            else
                sectionBreakColors(end + 1) = 2;
            end
        end%iTrial
    end%if

    % Plot lines dividing sections if necessary
    yL = ylim();
    if ~isempty(sectionBreakTrials)
        for iColor = 1:size(sectionBreakColorList, 1)
            currColor = sectionBreakColorList(iColor, :);
            currSectionBreaks = sectionBreakTrials(sectionBreakColors(iColor));
            for iSection = 1:numel(currSectionBreaks)
                currTrialNum = currSectionBreaks(iSection);                
                xx = [currTrialNum currTrialNum] - 0.5;
                yy = [0 ceil(max(tuningAmp(:)))];
                if iColor == 1 && iSection == 1
                    plotLines = plot(xx, yy, 'color', currColor, 'linewidth', 1.5);
                elseif iSection == 1
                    plotLines(end + 1) = plot(xx, yy, 'color', currColor, 'linewidth', 1.5);
                else
                    plot(xx, yy, 'color', currColor, 'linewidth', 1.5);
                end
            end
        end
        if iPlot == 1
            lgd = legend(plotLines, sectionColorLabels, 'autoupdate', 'off', 'location', 'none')
            lgd.FontSize = 10;
            lgd.Orientation = 'horizontal';
            lgd.Box = 'off';
            drawnow()
%             lgd.Position(3) = lgd.Position(3) * 1.2;
            firstAx = ax;
        end
    end

    if iPlot == nPlots
        xlabel('Trial number', 'fontsize', 13)
    else
        ax.XTickLabel = [];
    end
    xL = xlim;
    xlim([0, size(tuningAmp, 2) + 1])
    ylabel(num2str(iPlot))
end

for iAx = 1:numel(f.Children)
   if strcmp(f.Children(iAx).Tag, 'subaxis')
       f.Children(iAx).YLim = [0, 1.2 * max(tuningAmp(:))];
   end
end

% Plot title at top of figure
figTitleStr = {[bl.expID, ' - EB wedge visual tuning'], ...
        [figTitleText, '  (max - min)']};
if omitMoveVols
    figTitleStr{2} = [figTitleStr{2}, '  -  excl. move.'];
end
h = suptitle(figTitleStr);
h.FontSize = 12;

axPos = firstAx.Position;
lgd.Position(1) = axPos(1);
lgd.Position(2) = sum(axPos([2 4])) * 1.002;

% Save figure
if saveFig
   saveDir = fullfile(analysisDir, bl.expID);
   save_figure(f, saveDir, saveFileName);
end

catch ME; rethrow(ME); end

%% =================================================================================================
%% Plot summary of tuning curve amplitudes
%===================================================================================================

saveFig = 0;

omitMoveVols = 1; 

smWin = 3;

sourceData = bl.wedgeRawFlArr;
sourceData = bl.wedgeDffArr;
% sourceData = bl.wedgeZscoreArr;
sourceData = bl.wedgeExpDffArr;

try
% Generate figure labels and save file name
if isequal(sourceData, bl.wedgeDffArr)
    figTitleText = 'dF/F';
    saveFileName = 'tuning_curve_amplitude_summary_dff';
elseif isequal(sourceData, bl.wedgeRawFlArr)
    figTitleText = 'Raw F';
    saveFileName = 'tuning_curve_amplitude_summary_rawF';
elseif isequal(sourceData, bl.wedgeZscoreArr)
    figTitleText = 'Z-scored dF/F';
    saveFileName = 'tuning_curve_amplitude_summary_zscored-dff';
elseif isequal(sourceData, bl.wedgeExpDffArr)
    figTitleText = 'full exp dF/F';
    saveFileName = 'tuning_curve_amplitude_summary_exp-dff';
else
    warndlg('Warning: sourceData mismatch');
    figTitleText = '';
    saveFileName = 'tuning_curve_amplitude_summary_unknown';
end


% Replace volumes when fly was moving with nan if necessary
if omitMoveVols
   moveDistThreshVols = round(moveDistThresh * bl.volumeRate);
   moveDistArr = permute(repmat(bl.moveDistVols, 1, 1, 8), [1 3 2]);
   sourceData(moveDistArr < moveDistThreshVols) = nan; 
   saveFileName = [saveFileName, '_no-movement'];
end

% Get mean panels pos data
panelsPosVols = [];
for iVol = 1:size(bl.wedgeDffArr, 1)
    [~, currVol] = min(abs(bl.panelsFrameTimes - bl.volTimes(iVol)), [], 'omitnan');
    panelsPosVols(iVol) = bl.panelsPosX(currVol);
end
tuningData = [];
for iPos = 1:numel(unique(bl.panelsPosX))
    tuningData(iPos, :, :) = ...
            mean(sourceData(panelsPosVols == (iPos - 1), :, :), 1, 'omitnan'); % --> [barPos, wedge, trial]    
end

% Get mean panels pos data
panelsPosVols = [];
for iVol = 1:size(bl.wedgeDffArr, 1)
    [~, currVol] = min(abs(bl.panelsFrameTimes - bl.volTimes(iVol)));
    panelsPosVols(iVol) = bl.panelsPosX(currVol);
end
tuningData = [];
for iPos = 1:numel(unique(bl.panelsPosX))
    tuningData(iPos, :, :) = ...
            mean(sourceData(panelsPosVols == (iPos - 1), :, :), 1, 'omitnan'); % --> [barPos, wedge, trial]    
end

% Smooth data, then find the max and min for each trial and wedge
tuningDataSm = repeat_smooth(tuningData, 10, 'smWin', smWin);
minVals = []; maxVals = [];
for iTrial = 1:size(tuningDataSm, 3)
    minVals(:, iTrial) = min(tuningDataSm(:, :, iTrial), [], 1, 'omitnan'); % --> [wedge, trial]
    maxVals(:, iTrial) = max(tuningDataSm(:, :, iTrial), [], 1, 'omitnan'); % --> [wedge, trial]
end
xTickLabel = 1:size(minVals, 2);

% Cut out opto stim trials
optoStimTrials = [bl.usingOptoStim];
postStimTrials = [0, optoStimTrials(1:end - 1)];
postStimTrials = postStimTrials(~optoStimTrials);
optoStimPanels = bl.usingPanels;
prevStimPanels = [0, optoStimPanels(1:end - 1)];
prevStimPanels = prevStimPanels(~optoStimTrials);
xTickLabel = xTickLabel(~optoStimTrials);
minVals = minVals(:, ~optoStimTrials);  % --> [wedge, trial]
maxVals = maxVals(:, ~optoStimTrials);  % --> [wedge, trial]
tuningAmp = maxVals - minVals;          % --> [wedge, trial]


% Group by experiment epoch
blockAmps = [];
blockStartTrials = [1, find(postStimTrials)];
for iBlock = 1:numel(blockStartTrials)
    startTrial = blockStartTrials(iBlock);
    if iBlock == numel(blockStartTrials)
        endTrial = size(tuningAmp, 2);
    else
        endTrial = blockStartTrials(iBlock + 1) - 1;
    end
    blockAmps(:, iBlock) = mean(tuningAmp(:, startTrial:endTrial), 2); % --> [wedge, block]
end

% Pull out the tuning curve amplitudes for trials before and after opto stims
blockAmps = [];
blockStartTrials = find(postStimTrials);
for iTrial = 1:numel(blockStartTrials)
   preStim = blockStartTrials(iTrial) - 1;
   postStim = blockStartTrials(iTrial);
   blockAmps(:, :, iTrial) = tuningAmp(:, [preStim, postStim]); % --> [wedge, trial, stim trial num]
end
blockAmps = permute(blockAmps, [2, 1 3]); % --> [trial, wedge, stim trial Num]

% Create figure
nStimTrials = size(blockAmps, 3);
f = figure(1);clf;
f.Color = [1 1 1];
f.Position = [f.Position(1:2), nStimTrials * 300, 450];
usingPanels = bl.usingPanels(logical(bl.usingOptoStim));
optoStimTrials = bl.trialNum(logical(bl.usingOptoStim));
for iTrial = 1:nStimTrials
    subaxis(1, nStimTrials, iTrial, 'mr', 0.04, 'ml', 72 / f.Position(3), 'mt', 0.15, 'mb', 0.11); hold on; box off
    plot(blockAmps(:, :, iTrial), 'o-')
    ax = gca;
    ax.XTick = [1 2];
    ax.XTickLabel = {'Before', 'After'};
    ax.FontSize = 12;
    xlim([0.75 2.25]);
    if usingPanels(iTrial)
        title(['Opto + bar (trial ', num2str(optoStimTrials(iTrial)), ')'], 'FontSize', 12)
    else
        title(['Opto only (trial ', num2str(optoStimTrials(iTrial)), ')'], 'FontSize', 12)
    end
    
    % Plot mean amplitude for each group
    plot([mean(blockAmps(1, :, iTrial)), mean(blockAmps(2, :, iTrial))], ...
            'o-', 'Color', 'k', 'linewidth', 5)
    
    % Track Y limits of each axis
    if iTrial == 1
       figYLims = ylim();
       ylabel('Max - min tuning curve dF/F', 'fontsize', 13) 
    else
       currYLims = ylim();
       figYLims = [min([figYLims(1), currYLims(1)]), max([figYLims(2), currYLims(2)])];
    end
end

% Make each plot use the same Y limits
for iAx = 1:numel(f.Children)
    try
        f.Children(iAx).YLim = figYLims;
    catch
    end
end

figTitleStr = ([bl.expID, '  -   ', figTitleText]);
if omitMoveVols
    figTitleStr = [figTitleStr, '  -  excl. vol. within ', num2str(moveDistThresh), ...
            ' sec of movement'];
end
suptitle(figTitleStr)


% Save figure
if saveFig
   saveDir = fullfile(analysisDir, bl.expID);
   save_figure(f, saveDir, saveFileName);
end

catch ME; rethrow(ME); end

%% =================================================================================================
%% Plot average tuning curves for entire experiment
%===================================================================================================

saveFig = 0;

omitMoveVols = 1;

subtractMean = 0;

smWin = 3;


% sourceData = bl.wedgeRawFlArr;
% sourceData = bl.wedgeDffArr;
% sourceData = bl.wedgeZscoreArr;
sourceData = bl.wedgeExpDffArr;

% sourceData =  bl.dffArr;
% sourceData =  bl.rawFlArr;
% sourceData =  bl.zscoreArr;
% sourceData = bl.expDffArr;

try
roiNames = {mD(1).roiData.name};%{'Top', 'L1', 'L2', 'L3', 'Bottom', 'R1', 'R2', 'R3', 'Mean'};

%     sourceData = permute(sourceData(:, 1:end-1, :), [1 3 2]);
%     goodVolsRep = repmat(goodVols, 1, 1, size(sourceData, 3));
%     sourceData(~goodVolsRep) = nan;
%     sourceData = permute(sourceData, [1 3 2]);
%     roiNames(end) = [];


figSize = [];
figSize = [1100 850];

% Generate figure labels and save file name
if isequal(sourceData, bl.wedgeDffArr) || isequal(sourceData, bl.dffArr)
    figTitleText = 'dF/F';
    saveFileName = 'avg_tuning_curves_dff';
    yLabelText = 'dF/F';
elseif isequal(sourceData, bl.wedgeRawFlArr) || isequal(sourceData, bl.rawFlArr)
    figTitleText = 'raw F';
    saveFileName = 'avg_tuning_curves_rawF';
    yLabelText = 'Raw F';
elseif isequal(sourceData, bl.wedgeZscoreArr) || isequal(sourceData, bl.zscoreArr)
    figTitleText = 'Z-scored dF/F';
    saveFileName = 'avg_tuning_curve_zscored-dff';
    yLabelText = 'Z-scored dF/F';
elseif isequal(sourceData, bl.wedgeExpDffArr) || isequal(sourceData, bl.expDffArr)
    figTitleText = 'full exp dF/F';
    saveFileName = 'avg_tuning_curves_exp-dff';
    yLabelText = 'dF/F';
else
%     warndlg('Warning: sourceData mismatch');
    figTitleText = '';
    saveFileName = 'avg_tuning_curves_unknown';
end

%     figTitleText = 'Z-scored dF/F';
%     saveFileName = 'avg_tuning_curve_zscored-dff';
%     yLabelText = 'Z-scored dF/F';

% Replace volumes when fly was moving with nan if necessary
if omitMoveVols
   moveDistThreshVols = round(moveDistThresh * bl.volumeRate);
   moveDistArr = permute(repmat(bl.moveDistVols, 1, 1, size(sourceData, 2)), [1 3 2]);
   sourceData(moveDistArr < moveDistThreshVols) = nan; 
   saveFileName = [saveFileName, '_no-movement'];
end

% Get mean panels pos data
panelsPosVols = [];
for iVol = 1:size(sourceData, 1)
    [~, currVol] = min(abs(bl.panelsFrameTimes - bl.volTimes(iVol)));
    panelsPosVols(iVol) = bl.panelsPosX(currVol);
end
plotData = [];
for iPos = 1:numel(unique(bl.panelsPosX))
    plotData(iPos, :, :) = ...
            mean(sourceData(panelsPosVols == (iPos - 1), :, :), 1, 'omitnan'); % --> [barPos, wedge, trial]    
end

% Shift data so center of plot is directly in front of the fly
shiftData = cat(1, plotData(92:96, :, :), plotData(1:91, :, :));
shiftDataSm = repeat_smooth(shiftData, 10, 'smWin', smWin);

% Offset data so that each plot is centered at zero
shiftDataOffset = [];
for iTrial = 1:size(shiftDataSm, 3)
   for iWedge = 1:size(shiftDataSm, 2)
      currData = shiftDataSm(:, iWedge, iTrial); % --> [barPos]
      shiftDataOffset(:, iWedge, iTrial) = currData - mean(currData, 'omitnan');
   end    
end

% Find max and min values 
% yMax = max(shiftDataOffset(800:end));
% yMin = min(shiftDataOffset(800:end));
yMax = max(shiftDataOffset(:), [], 'omitnan');
yMin = min(shiftDataOffset(:), [], 'omitnan');
range = max(abs([yMax, yMin]), [], 'omitnan') *  2.5;

% Create figure and plots
f = figure(11);clf;
f.Color = [1 1 1];
if ~isempty(figSize)
   f.Position(3:4) = figSize; 
end
hold on;
nROIs = size(shiftDataOffset, 2);
cm = jet(nROIs);
cm(end + 1, :) = 0;
allROIMean = squeeze(mean(mean(shiftDataOffset, 3, 'omitnan'), 2, 'omitnan'));
if subtractMean
    n = nROIs;
else
    n = nROIs + 1;
    roiNames{end + 1} = 'Mean';
end
for iROI = 1:n
    
    % Plot data
    plotX = -180:3.75:(180 - 3.75);
    if iROI <= nROIs
        if subtractMean
            plotData = squeeze(mean(shiftDataOffset(:, iROI, :), 3, 'omitnan')) - allROIMean;
        else
            plotData = squeeze(mean(shiftDataOffset(:, iROI, :), 3, 'omitnan'));
        end
    else
        plotData = allROIMean;
    end
    plot(plotX, plotData, 'linewidth', 2, 'color', cm(iROI, :));
    
end
legend(roiNames, 'location', 'eastoutside');
xlabel('Bar position (degrees from front of fly)', 'fontsize', 14);
ylabel(yLabelText, 'fontsize', 14);


% Add figure title
figTitleStr = {[bl.expID, '  -  Visual tuning (mean ', figTitleText, ')']};
if subtractMean
   figTitleStr{2} = 'Mean tuning curve of all ROIs subtracted'; 
end
if omitMoveVols
    figTitleStr{1} = [figTitleStr{1}, '  -  excluding volumes within ', num2str(moveDistThresh), ...
            ' sec of movement'];
end
h = suptitle(figTitleStr);
h.FontSize = 14;


% Save figure
if saveFig
    if subtractMean
        saveFileName = [saveFileName, '_meanSub'];
    end
    saveDir = fullfile(analysisDir, bl.expID);
    save_figure(f, saveDir, saveFileName);
end

catch ME; rethrow(ME); end






%% =================================================================================================
%% BUMP AMPLITUDE COMPARISONS
%===================================================================================================

%% Within-wedge amplitude for each trial 

iTrial = 5;
smWin = 5;

% for iTrial = 1:size(test, 3)
currTrialData = bl.wedgeExpDffArr(:, :, iTrial);
smData = smoothdata(currTrialData, 1, 'gaussian', smWin);


% Identify 5th and 95th percentile for each wedge
smDataSorted = sort(smData);
pct5vals = smDataSorted(round(size(smDataSorted, 1) * 0.05), :);
pct95vals = smDataSorted(round(size(smDataSorted, 1) * 0.95), :);

bumpAmp = pct95vals - pct5vals;

% end

f = figure(1); clf;
ax = subaxis(1, 2, 1);
imagesc(smData');
ax.YTickLabel = round(bumpAmp, 2);

subaxis(1, 2, 2)
plot(smDataSorted);


%% Instantaneous bump amplitude across wedges for each time point 

smWin = 5;
saveFig = 0;

smData = smoothdata(bl.wedgeExpDffArr, 1, 'gaussian', smWin);
smDataExpSort = sort(reshape(permute(smData, [2 1 3]), size(smData, 2), []), 2);
wedgeDff95pct = smDataExpSort(:, round(size(smDataExpSort, 2) * 0.95));
smDataNorm = smData ./ repmat(wedgeDff95pct', size(smData, 1), 1, size(smData, 3));

plotData = smData;
plotData = smDataNorm;

bumpAmp = [];
for iTrial = 1:size(bl.wedgeExpDffArr, 3)
    currTrialData = plotData(:, :, iTrial);
    bumpAmp(:, iTrial) = max(currTrialData, [], 2) - min(currTrialData, [], 2);
end


% Plot mean bump amplitude for each entire trial compared to overall mean dF/F
f = figure(2); clf
subaxis(2, 1, 1)
plot(mean(bumpAmp, 1), '-o');
hold on;
plot(squeeze(mean(max(plotData, [], 2), 1)), '-o')
plot(squeeze(mean(min(plotData, [], 2), 1)), '-o')
legend({'Mean amp', 'Avg max', 'Avg min'})

title('Amplitude');

subaxis(2, 1, 2)
plot(multi_mean(plotData, [1 2]), '-o');
title('Mean dF/F')

%
figure(9);clf;
imagesc(reshape(smData, size(smData, 1), [])');
figure(10);clf;
imagesc(reshape(smDataNorm, size(smData, 1), [])');

% ------------------------------------------------------------------------------

% Plot volume-by-volume amplitude throughout the entire experiment
f = figure(8);clf; set(gcf, 'color', [1 1 1])
ax(1) = subaxis(16, 1, 1:9, 'mr', 0.04, 'ml', 0.05, 'sv', 0.08);
hold on
xx = (1:numel(bumpAmp)) / bl.volumeRate;
plot(xx, bumpAmp(:), 'color', [0 0 0.9])
xlim([xx(1) xx(end)])
ax(1).XTickLabel = 1:size(smData, 3);
ax(1).FontSize = 12;
ax(1).YLabel.String = 'Bump amplitude';

% plot(xx, repeat_smooth(bumpAmp(:), 500, 'smWin', 10), 'color', 'r', 'linewidth', 1.5)

% Add lines to mark trial boundaries and text labels for trial number
yL = ylim();
xTickLocs = [];
for iTrial = 2:size(smData, 3)
    xVal = bl.trialDuration * (iTrial - 1);
    plot([xVal, xVal], [0 yL(2)], ':', 'color', 'k', 'linewidth', 1);
    xTickLocs(end + 1) = xVal;
end
xTickLocs(end + 1) = xTickLocs(end) + bl.trialDuration;
ax(1).XTick = xTickLocs - (bl.trialDuration / 2);
ax(1).TickLength = [0 0];

% Add annotation lines
expSectionBreaks = ((sectionBreakTrials - 1) * bl.trialDuration);
 if ~isempty(expSectionBreaks)
        for iColor = 1:size(sectionBreakColorList, 1)
            currColor = sectionBreakColorList(iColor, :);
            currSectionBreaks = expSectionBreaks(sectionBreakColors(iColor));
            for iSection = 1:numel(currSectionBreaks)
                currTrialNum = currSectionBreaks(iSection);                
                xx = [currTrialNum currTrialNum] - (0 * bl.trialDuration);
                yy = [0 yL(2)];
                if iColor == 1 && iSection == 1
                    plotLines = plot(xx, yy, 'color', currColor, 'linewidth', 3);
                elseif iSection == 1
                    plotLines(end + 1) = plot(xx, yy, 'color', currColor, 'linewidth', 3);
                else
                    plot(xx, yy, 'color', currColor, 'linewidth',3);
                end
            end
        end
        

        lgd = legend(plotLines, sectionColorLabels(1:numel(plotLines)), 'autoupdate', 'off', ...
                'location', 'northoutside');
        lgd.FontSize = 10;
        lgd.Orientation = 'horizontal';
        lgd.Box = 'off';
        drawnow()
 end
ylim(yL);


% Plot the full data underneath
ax(2) = subaxis(16, 1, 10:16);
xx = (1:numel(bumpAmp)) / bl.volumeRate;
imagesc([0 xx(end)], [0.5, size(plotData, 2) + 0.5], reshape(permute(plotData, [2 1 3]), ...
        size(plotData, 2), []));
linkaxes(ax, 'x');
ax(2).YTickLabel = [];
ax(2).FontSize = 12;
ax(2).YLabel.String = 'Normalized dF/F';
ax(2).XLabel.String = 'Time (s)';

suptitle(bl.expID)

if saveFig
   saveDir = fullfile(analysisDir, bl.expID);
   save_figure(f, saveDir, 'bump_amplitude_normalized');
end







