
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\20191205-2_38A11_ChR_60D05_7f\ProcessedData';
analysisDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\Analysis';
load(fullfile(parentDir, 'analysis_data.mat'));

%% COMBINE DATA FROM A BLOCK OF COMPATIBLE TRIALS

blTrials = [1:24];

bD = mD(ismember([mD.trialNum], blTrials));

% Check that have compatible values in scalar fields 
if numel(unique([bD.daqSampRate])) > 1 || ...
        numel(unique({bD.expID})) > 1 || ...
        numel(unique([bD.nDaqSamples])) > 1 || ...
        numel(unique([bD.nPlanes])) > 1 || ...
        numel(unique([bD.nPanelsFrames])) > 1 || ...
        numel(unique([bD.nVolumes])) > 1 || ...
        numel(unique([bD.panelsDisplayRate])) > 1 || ...
        numel(unique([bD.trialDuration])) > 1 || ...
        numel(unique([bD.using2P])) > 1 || ...
        numel(unique([bD.volumeRate])) > 1 || ...
        numel(unique([bD.nVolumes])) > 1 || ...
        numel(unique([bD.panelsCycleFrames])) > 1 || ...
        numel(unique([bD.panelsCycleTime])) > 1
    disp('Error: the specified trials are not compatible');
end

% Add data that applies to all trials
bl = [];
bl.daqSampRate = bD(1).daqSampRate;
bl.daqSampTimes = bD(1).daqSampTimes;
bl.expID = bD(1).expID;
bl.nDaqSamples = bD(1).nDaqSamples;
bl.nPanelsFrames = bD(1).nPanelsFrames;
bl.nPlanes = bD(1).nPlanes;
bl.nVolumes = bD(1).nVolumes;
bl.optoStimTiming = bD(1).optoStimTiming;
bl.panelsCycleFrames = bD(1).panelsCycleFrames;
bl.panelsCycleTime = bD(1).panelsCycleTime;
bl.panelsDisplayRate = bD(1).panelsDisplayRate;
bl.panelsFrameTimes = bD(1).panelsFrameTimes;
bl.panelsPattern = bD(1).panelsPattern;
bl.panelsPosX = bD(1).panelsPosX;
bl.panelsPosY = bD(1).panelsPosY;
bl.refImages = bD(1).refImages;
bl.trialDuration = bD(1).trialDuration;
bl.volTimes = bD(1).volTimes;


% Add compiled data
bl.dffArr = []; bl.dffVectAvgRad = []; bl.dffVectAvgWedge = []; bl.dffVectStrength = [];
bl.ftData = []; bl.optoStimOnsetTimes = []; bl.optoStimOffsetTimes = []; bl.rawFlArr = [];
bl.wedgeDffArr = []; bl.wedgeZscoreArr = []; bl.zscoreArr = []; bl.usingOptoStim = [];
bl.trialNum = []; bl.usingPanels = []; bl.wedgeRawFlArr = [];
for iTrial = 1:numel(bD)
    bl.dffArr(:, :, iTrial) = bD(iTrial).dffMat;                        % --> [volume, glom, trial]
    bl.rawFlArr(:, :, iTrial) = bD(iTrial).rawFlMat;                    % --> [volume, glom, trial]
    bl.dffVectAvgRad(:, iTrial) = bD(iTrial).dffVectAvgRad;             % --> [volume, trial]
    bl.dffVectAvgWedge(:, iTrial) = bD(iTrial).dffVectAvgWedge;         % --> [volume, trial]        
    bl.dffVectStrength(:, iTrial) = bD(iTrial).dffVectStrength;         % --> [volume, trial]
    bl.optoStimOnsetTimes{iTrial} = bD(iTrial).optoStimOnsetTimes;      % --> {trial}
    bl.optoStimOffsetTimes{iTrial} = bD(iTrial).optoStimOffsetTimes;    % --> {trial}
    bl.usingOptoStim(iTrial) = bD(iTrial).usingOptoStim;                % --> [trial]
    bl.wedgeRawFlArr(:, :, iTrial) = bD(iTrial).wedgeRawFlMat;          % --> [volume, wedge, trial]
    bl.wedgeDffArr(:, :, iTrial) = bD(iTrial).wedgeDffMat;              % --> [volume, wedge, trial]
    bl.wedgeZscoreArr(:, :, iTrial) = bD(iTrial).wedgeZscoreMat;        % --> [volume, wedge, trial]
    bl.zscoreArr(:, :, iTrial) = bD(iTrial).zscoreMat;                  % --> [volume, glom, trial]
    bl.trialNum(iTrial) = bD(iTrial).trialNum;                          % --> [trial]
    bl.usingPanels(iTrial) = bD(iTrial).usingPanels;                    % --> [trial]
    
    % Trim FicTrac frames if necessary 
    targetFrames = min(cellfun(@numel, {bD.ftFrameTimes}));
    if numel(bD(iTrial).ftFrameTimes) == targetFrames
       bl.ftFrameTimes = bD(iTrial).ftFrameTimes; 
    end
    % All are --> [frame, trial]
    bl.ftData.intX(:, iTrial) = bD(iTrial).ftData.intX(1:targetFrames);
    bl.ftData.intY(:, iTrial) = bD(iTrial).ftData.intY(1:targetFrames);
    bl.ftData.intHD(:, iTrial) = bD(iTrial).ftData.intHD(1:targetFrames);
    bl.ftData.moveSpeed(:, iTrial) = bD(iTrial).ftData.moveSpeed(1:targetFrames);
    bl.ftData.intFwMove(:, iTrial) = bD(iTrial).ftData.intFwMove(1:targetFrames);
    bl.ftData.intSideMove(:, iTrial) = bD(iTrial).ftData.intSideMove(1:targetFrames);
    bl.ftData.yawSpeed(:, iTrial) = bD(iTrial).ftData.yawSpeed(1:targetFrames);
    bl.ftData.fwSpeed(:, iTrial) = bD(iTrial).ftData.fwSpeed(1:targetFrames);
    bl.ftData.sideSpeed(:, iTrial) = bD(iTrial).ftData.sideSpeed(1:targetFrames);
    
    bl.ftData.moveSpeed(1,iTrial) = bl.ftData.moveSpeed(2,iTrial);
end

% Sort fields alphabetically
bl = orderfields(bl);


%% ===================================================================================================
%% Plot tuning curves for each wedge across trials
%===================================================================================================


saveFig = 1;

sourceData = bl.wedgeRawFlArr;
% sourceData = bl.wedgeDffArr;
% sourceData = bl.wedgeZscoreArr;

showStimTrials = 1;



% Generate figure labels and save file name
if isequal(sourceData, bl.wedgeDffArr)
    figTitleText = 'dF/F';
    saveFileName = '2D_tuning_curves_dff';
elseif isequal(sourceData, bl.wedgeRawFlArr)
    figTitleText = 'raw F';
    saveFileName = '2D_tuning_curves_rawF';
elseif isequal(sourceData, bl.wedgeZscoreArr)
    figTitleText = 'Z-scored dF/F';
    saveFileName = '2D_tuning_curves_zscored-dff';
else
    errordlg('Error: sourceData mismatch');
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
            mean(sourceData(panelsPosVols == (iPos - 1), :, :), 1); % --> [barPos, wedge, trial]    
end
meanDataSmooth = smoothdata(meanData, 1, 'gaussian', 3);

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
            optoStimTrials = find([mD.usingOptoStim]);
            for iTrial = 1:numel(optoStimTrials)
                plotY = optoStimTrials(iTrial) - skippedTrials - 0.5;
                skippedTrials = skippedTrials + 1;
                if mD(optoStimTrials(iTrial)).usingPanels
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
h = suptitle({[bl.expID, '  -  Visual tuning of EB wedges (mean ', figTitleText, ')'], ...
        'Green line  =  trial with opto + visual stim', ...
        '  Red line  = trial with opto stim only'});
h.FontSize = 14;
    
% Save figure
if saveFig
   saveDir = fullfile(analysisDir, bl.expID);
   save_figure(f, saveDir, saveFileName);
end

%% ===================================================================================================    
%% Plot as lines instead of using imagesc
% ===================================================================================================

saveFig = 1;

smWin = 2;
sourceData = bl.wedgeRawFlArr;
% sourceData = bl.wedgeDffArr;
% sourceData = bl.wedgeZscoreArr;

figSize = [];
figSize = [1250 980];

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
else
    errordlg('Error: sourceData mismatch');
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
            mean(sourceData(panelsPosVols == (iPos - 1), :, :), 1); % --> [barPos, wedge, trial]    
end

% Shift data so center of plot is directly in front of the fly
shiftData = cat(1, plotData(92:96, :, :), plotData(1:91, :, :));
shiftDataSm = repeat_smooth(shiftData, 10, 'smWin', smWin);

% Offset data so that each plot is centered at zero
shiftDataOffset = [];
for iTrial = 1:size(shiftDataSm, 3)
   for iWedge = 1:size(shiftDataSm, 2)
      currData = shiftDataSm(:, iWedge, iTrial); % --> [barPos]
      shiftDataOffset(:, iWedge, iTrial) = currData - mean(currData);
   end    
end

% Find max and min values 
% yMax = max(shiftDataOffset(800:end));
% yMin = min(shiftDataOffset(800:end));
yMax = max(shiftDataOffset(:));
yMin = min(shiftDataOffset(:));
range = max(abs([yMax, yMin])) *  2.5;

% Create figure and plots
f = figure(1);clf;
f.Color = [1 1 1];
if ~isempty(figSize)
   f.Position(3:4) = figSize; 
end
for iPlot = 1:nPlots
    ax = subaxis(1, nPlots, iPlot, 'mt', 0, 'ml', 0.05, 'mr', 0.03, 'mb', 0.08, 'sh', 0.05);
    hold on;
    currData = squeeze(shiftDataOffset(:, iPlot, :)); % --> [barPos, trial]
    currData = currData(:, ~[bl.usingOptoStim]);
    
    % Plot data
    plotX = -180:3.75:(180 - 3.75);
    
    % Separate data into distinct rows
    offsets = [];
    currData = fliplr(currData); % Flip so first trial is in the last column
    for iTrial = 1:size(currData, 2)
        currOffset = (range * (iTrial - 1));
        offsets(iTrial) = currOffset;
        currData(:, iTrial) = currData(:, iTrial) + currOffset;
        plot([plotX(1), plotX(end)], [currOffset, currOffset], '--', 'color', 'b') % Plot zero lines
    end
    plot(plotX, currData, 'color', 'k', 'linewidth', 1.25);
    yL(1) = yMin * 1.25;
    yL(2) = range * (size(currData, 2));
    hold on; plot([0 0], yL, 'color', 'b', 'linewidth', 1)
    ylim(yL);
    
    xlabel('Bar position (deg)');
    if iPlot == 1
        ylabel('Trial', 'fontsize', 16);
    end
    ax.YTick = offsets;
    yTickLabel = size(shiftDataOffset, 3):-1:1;
    ax.YTickLabel = yTickLabel(~bl.usingOptoStim);
    
    % Indicate missing opto stim trials
    skippedTrials = 0;
    optoStimTrials = find([bl.usingOptoStim]);   
    for iTrial = 1:numel(optoStimTrials)
        offsetsRev = offsets(end:-1:1);
        plotOffset = offsetsRev(optoStimTrials(iTrial) - skippedTrials - 1) - (range / 2);
        skippedTrials = skippedTrials + 1;
        if bl.usingPanels(optoStimTrials(iTrial))
            lineColor = 'green';
        else
            lineColor = 'red';
        end
        plot([plotX(1), plotX(end)], [plotOffset, plotOffset], '-', 'color', lineColor, ...
                'linewidth', 1.5)
    end
    
end

% Add figure title
h = suptitle({[bl.expID, '  -  Visual tuning of EB wedges (mean ', figTitleText, ')'], ...
        'Green line  =  trial with opto + visual stim', ...
        '  Red line  = trial with opto stim only'});
h.FontSize = 14;
    
% Save figure
if saveFig
   saveDir = fullfile(analysisDir, bl.expID);
   save_figure(f, saveDir, saveFileName);
end

%% ===================================================================================================
%% Plot min and max values from the tuning curves
%===================================================================================================

saveFig = 1;

smWin = 2;

sourceData = bl.wedgeRawFlArr;
% sourceData = bl.wedgeDffArr;
% sourceData = bl.wedgeZscoreArr;
    
figSize = [];
figSize = [400 925];

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
else
    errordlg('Error: sourceData mismatch');
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
            mean(sourceData(panelsPosVols == (iPos - 1), :, :), 1); % --> [barPos, wedge, trial]    
end

% Smooth data, then find the max and min for each trial and wedge
tuningDataSm = repeat_smooth(tuningData, 10, 'smWin', smWin);
minVals = []; maxVals = [];
for iTrial = 1:size(tuningDataSm, 3)
    minVals(:, iTrial) = min(tuningDataSm(:, :, iTrial), [], 1); % --> [wedge, trial]
    maxVals(:, iTrial) = max(tuningDataSm(:, :, iTrial), [], 1); % --> [wedge, trial]
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
   f.Position(3:4) = figSize; 
end
globalMax = 0;
for iPlot = 1:nPlots
    ax = subaxis(nPlots, 1, iPlot, 'mt', 0.06, 'mb', 0.06, 'sv', 0.03, 'mr', 0.03, 'ml', 0.08);
    hold on;
%     plot(1:size(minVals, 2), minVals(iPlot, :), 'o', 'color', 'b')
%     plot(1:size(maxVals, 2), maxVals(iPlot, :), 'o', 'color', 'm')
    plot(1:size(tuningAmp, 2), tuningAmp(iPlot, :), 'o', 'color', 'k')
    stimTrials = find(postStimTrials);
    for iTrial = stimTrials
        if prevStimPanels(iTrial)
            lineColor = 'green';
        else
            lineColor = 'red';
        end
        yL = ylim();
        plot([iTrial - 0.5, iTrial - 0.5], yL, '-', 'color', lineColor, ...
            'linewidth', 1.5)
        globalMax = max([globalMax, yL(2)]);
    end
    if iPlot == nPlots
        xlabel('Trial number', 'fontsize', 13)
    else
        ax.XTickLabel = [];
    end
    xL = xlim;
    xlim([0, xL(2) + 1])
%     ylabel(num2str(iPlot))
end

% for iAx = 1:numel(f.Children)
%    if strcmp(f.Children(iAx).Tag, 'subaxis')
%        f.Children(iAx).YLim = [0 globalMax]
%    end
% end

% Plot title at top of figure
h = suptitle({[bl.expID, ' - EB wedge visual tuning'], ...
        [figTitleText, '  (max - min)'], ... 
        'Green line  =  trial with opto + visual stim', ...
        '  Red line  = trial with opto stim only'});
h.FontSize = 12;

% Save figure
if saveFig
   saveDir = fullfile(analysisDir, bl.expID);
   save_figure(f, saveDir, saveFileName);
end

%% ===================================================================================================
%% Plot behavior as a function of bar position
%===================================================================================================

smReps = 6;
smWin = 5;


% Get mean panels pos data
panelsPosVols = [];
for iFrame = 1:numel(bl.ftFrameTimes)
    [~, currFrame] = min(abs(bl.panelsFrameTimes - bl.ftFrameTimes(iFrame)));
    panelsPosFrames(currFrame) = bl.panelsPosX(currFrame);
end
meanSpeed = []; meanYawSpeed = []; meanYawVel = [];
for iPos = 1:numel(unique(bl.panelsPosX))
    meanSpeed(iPos, :) = ...
            mean(bl.ftData.moveSpeed(panelsPosFrames == (iPos - 1), :), 1);    % --> [barPos, trial]    
    meanYawVel(iPos, :) = ...
            mean(bl.ftData.yawSpeed(panelsPosFrames == (iPos - 1), :), 1);     % --> [barPos, trial] 
    meanYawSpeed(iPos, :) = ...
            mean(abs(bl.ftData.yawSpeed(panelsPosFrames == (iPos - 1), :)), 1);% --> [barPos, trial]
end


% Shift to make center of plot directly in front of fly
meanSpeed = [meanSpeed(92:96, :); meanSpeed(1:91, :)];
meanYawVel = [meanYawVel(92:96, :); meanYawVel(1:91, :)];
meanYawSpeed = [meanYawSpeed(92:96, :); meanYawSpeed(1:91, :)];

% plotVar = meanSpeed;
plotVar = meanYawSpeed;
plotVar(plotVar > 0.8) = 0.8;
plotX = -180:3.75:(180 - 3.75);

figure(1);clf;
plotData = plotVar(:, ~[bl.usingOptoStim]);
imagesc(plotX, [1 size(plotData, 2)], ...
    repeat_smooth(plotData, smReps, 'smWin', smWin)');
hold on; plot([0 0], [0 size(plotData, 2) + 1], 'color', 'k', 'linewidth', 2)
ylim([0.5 size(plotData, 2) + 0.5])
% colormap('bluewhitered')

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