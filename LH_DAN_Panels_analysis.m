
startDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';
parentDir = uigetdir(startDir, 'Select an experiment directory');

analysisDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\Analysis';

load(fullfile(parentDir, 'analysis_data.mat'));


load(fullfile(parentDir, 'flowMags.mat'));


%%


currTrial = 7;
td = mD([mD.trialNum] == currTrial);

plotROIs = [];

smWin = 5;
moveThresh = 0.08;
avgWin = [3 4];

flData = td.dffMat;
% flData = td.rawFlMat;
% flData = td.zscoreMat;
flData = td.expDffMat;

flowYLim = [0 0.6];

try

% PLOT DATA THROUGHOUT TRIAL
figure(1);clf;
clear allAx

subaxis(3, 3, 1, 'mb', 0.05, 'mt', 0.03)
subaxis(3, 3, 1:3)
allAx(1) = gca;


% Panels bar position
if td.usingPanels
    plot(td.panelsFrameTimes, td.panelsPosX, 'linewidth', 1.5, 'color', 'b');%colorbar
    barCenteredFrames = find(td.panelsPosX == 44);
end
if td.usingOptoStim
   hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
end
ylabel('Bar position')
% allAx(end).YTickLabel = [];
% allAx(end).XTickLabel = [];

subaxis(3, 3, 4:6)
allAx(end + 1) = gca;
% Mean dF/F data for each ROI
if isempty(plotROIs)
   plotROIs = 1:numel(td.roiData);    
end
% plotData = repeat_smooth(flData(:, plotROIs), 10, 'smwin', 3);
plotData = smoothdata(squeeze(flData(:, plotROIs)), 1, 'gaussian', smWin);
plot(repmat(td.volTimes', 1, size(plotData, 2)), plotData, 'linewidth', 1); 
hold on;
if td.usingOptoStim
   hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
end
ylim([min(plotData(:)) max(plotData(:))])
xlim([0 td.trialDuration])
roiNames = {td.roiData(plotROIs).name};
legend(roiNames)
ylabel('dF/F');
    
subaxis(3, 3, 7:9)
allAx(end + 1) = gca;
% 
% Optic flow data to identify flailing
flowThresh = 0.08;
currFlow = meanFlowMags{[mD.trialNum] == currTrial};
currFlow(end) = 0;
flowFrameDur = median(diff(td.ftFrameTimes));
flowFrameTimes = (1:1:numel(currFlow)) * flowFrameDur; % NOTE: this is only an approximation
plotData = repeat_smooth(currFlow, 20, 'dim', 2, 'smwin', 6);
plotData = plotData - min(plotData); 
plot(flowFrameTimes, plotData, 'color', 'k');
hold on;
if td.usingOptoStim
   hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
end
plot([td.ftFrameTimes(1), td.ftFrameTimes(end)], [moveThresh, moveThresh],...
        'linewidth', 0.5, 'color', 'r');
ylim(flowYLim)
ylabel('Optic flow (flailing)')

% Link the X-axis limits across all plots
% linkaxes(allAx(3:end), 'x');
linkaxes(allAx(:), 'x');
xlim([0 td.trialDuration])


% Calculate average responses to each opto stim period
baselineDur = avgWin(1);
postStimDur = avgWin(2);

if td.usingOptoStim
   
    stimWinDffData = [];
    baselineVols = floor(baselineDur * td.volumeRate);
    postStimVols = floor(postStimDur * td.volumeRate);
    for iStim = 1:numel(td.optoStimOnsetTimes)
        
        [~, onsetVol] = min(abs(td.volTimes - td.optoStimOnsetTimes(iStim)));
        stimDurVols = floor((td.optoStimOffsetTimes(iStim) - td.optoStimOnsetTimes(iStim)) ...
                * td.volumeRate);
        
        % Separate data around stim times
        startVol = onsetVol - baselineVols;
        endVol = onsetVol + stimDurVols + postStimVols;
        if startVol > 0 && endVol < td.nVolumes
            stimWinDffData(:, :, iStim) = flData(startVol:endVol, :); % --> [volume, ROI, stim]
        end
        stimWinDffAvg = mean(stimWinDffData, 3); % --> [volume, wedge];
        
    end
    plotTimes = td.volTimes(1:size(stimWinDffAvg, 1)) - td.volTimes(baselineVols + 1);
    
    figure(2);clf;hold on;
    for iWedge = 1:numel(plotROIs)
        plot(plotTimes, smoothdata(stimWinDffAvg(:, plotROIs(iWedge)), 1, 'gaussian', smWin), ...
                'linewidth', 2);
%                 'color', cm(iWedge, :), 'linewidth', 2);
    end
    plot_stim_shading([0 (td.optoStimOffsetTimes(1) - td.optoStimOnsetTimes(1))])
    legend(roiNames)
end
catch foldME; rethrow(foldME); end


%% Identify epochs of quiescence vs. flailing

moveThresh = 0.08;

for iTrial = 1:numel(mD)

   currTrialFlow = meanFlowMags{iTrial};
   smFlow = repeat_smooth(currTrialFlow, 20, 'dim', 2, 'smwin', 6);
   smFlow = smFlow - min(smFlow(:));
   moveFrames = smFlow > moveThresh;
   moveFrames([1 numel(moveFrames)]) = 1; % Assume that fly was moving at start and end of trial
   
   % For each quiescence frame, find the distance to the nearest movement frame
   moveFrameDist = zeros(1, numel(moveFrames));
   moveFrameInds = find(moveFrames);
   for iFrame = 1:numel(moveFrames)
       moveFrameDist(iFrame) = min(abs(moveFrameInds - iFrame));
   end
   
   mD(iTrial).moveFlowFrames = moveFrames;
   mD(iTrial).moveFlowFrameDist = moveFrameDist;
   
end


%% Extract fluorescence data for visual stim epochs only when fly is quiescent

trialNums = [1:5];
panelsStimPositions = 14;%[19 67 43]; % 
stimPosNames = {'Full field flash'}; % {'Left bar', 'Right bar', 'Center bar'}; %
quiescenceWinSec = [4 6];

stimWinFlow = [];

try
stimWinFlData = struct();
for iTrial = 1:numel(trialNums)
    
    currTrialInd = find([mD.trialNum] == trialNums(iTrial));
    
    % Identify stim epochs for each panels stim position
    currBarPosData = mD(currTrialInd).panelsPosX;
    stimOnsetFrames = {}; stimOffsetFrames = {};
    for iPos = 1:numel(panelsStimPositions)
        stimOnsetFrames = [];
        stimOffsetFrames = [];
        if iTrial == 1
            stimWinFlData(iPos).rawFl = [];
            stimWinFlData(iPos).dff = [];
            stimWinFlData(iPos).expDff = [];
            stimWinFlData(iPos).zscore = [];
        end
        for iFrame = 2:(numel(currBarPosData)-1)
           if currBarPosData(iFrame - 1) ~= panelsStimPositions(iPos) && ...
                        currBarPosData(iFrame) == panelsStimPositions(iPos)
              stimOnsetFrames(end + 1) = iFrame;
           elseif currBarPosData(iFrame + 1) ~= panelsStimPositions(iPos) && ...
                        currBarPosData(iFrame) == panelsStimPositions(iPos)
              stimOffsetFrames(end + 1) = iFrame;      
           end
        end
        
        % Convert from panels frames to imaging volumes
        stimOnsetTimes = mD(currTrialInd).panelsFrameTimes(stimOnsetFrames);
        stimOffsetTimes = mD(currTrialInd).panelsFrameTimes(stimOffsetFrames);
        
        % Filter out stim epochs that are accompanied by movement
        currFlow = meanFlowMags{currTrialInd};
        currFlow(end) = 0;
        flowFrameDur = median(diff(mD(currTrialInd).ftFrameTimes));
        flowFrameTimes = (1:1:numel(currFlow)) * flowFrameDur; % NOTE: this is only an approximation
        
        stimOnsetFlowFrames = []; stimOffsetFlowFrames = [];
        for iStim = 1:numel(stimOnsetTimes)
           [~, stimOnsetFlowFrames(iStim)] = min(abs(stimOnsetTimes(iStim) - flowFrameTimes));
           [~, stimOffsetFlowFrames(iStim)] = min(abs(stimOffsetTimes(iStim) - flowFrameTimes));
        end
        quiescenceWinFlowFrames = ceil(quiescenceWinSec ./ flowFrameDur);
        stimWinStartFlowFrames = stimOnsetFlowFrames - quiescenceWinFlowFrames(1);
        stimWinEndFlowFrames = stimOffsetFlowFrames + quiescenceWinFlowFrames(2);        
        
        stimOnsetTimes(stimWinStartFlowFrames < 1) = [];
        stimOffsetTimes(stimWinStartFlowFrames < 1) = [];
        stimWinEndFlowFrames(stimWinStartFlowFrames < 1) = [];
        stimWinStartFlowFrames(stimWinStartFlowFrames < 1) = [];
        
        stimOnsetTimes(stimWinEndFlowFrames > ...
                numel(mD(currTrialInd).moveFlowFrameDist)) = [];
        stimOffsetTimes(stimWinEndFlowFrames > ...
                numel(mD(currTrialInd).moveFlowFrameDist)) = [];
        stimWinStartFlowFrames(stimWinEndFlowFrames > ...
                numel(mD(currTrialInd).moveFlowFrameDist)) = [];
        stimWinEndFlowFrames(stimWinEndFlowFrames > ...
                numel(mD(currTrialInd).moveFlowFrameDist)) = [];

            
        minMoveDists = [];
        for iStim = 1:numel(stimOnsetTimes)
            minMoveDists(iStim) = min(mD(currTrialInd).moveFlowFrameDist( ...
                    stimWinStartFlowFrames(iStim):stimWinEndFlowFrames(iStim)));
        end
        
        stimWinOnsetTimes = stimOnsetTimes(logical(minMoveDists)) - quiescenceWinSec(1);
        stimWinOffsetTimes = stimOffsetTimes(logical(minMoveDists)) + quiescenceWinSec(2);
         
        % Extract fluorescence data from target stim window for each ROI
        stimWinOnsetVols = floor(stimWinOnsetTimes .* mD(currTrialInd).volumeRate);
        stimWinOffsetVols = floor(stimWinOffsetTimes .* mD(currTrialInd).volumeRate);
        for iStim = 1:numel(stimWinOnsetVols)
            currTd = mD(currTrialInd);
            currStart = stimWinOnsetVols(iStim);
            currEnd = stimWinOffsetVols(iStim);
            if ~isempty(stimWinFlData(iPos).rawFl) && ...
                        (currEnd - currStart + 1) > size(stimWinFlData(iPos).rawFl, 1)
               currEnd = currEnd - 1;
            elseif ~isempty(stimWinFlData(iPos).rawFl) && ...
                        (currEnd - currStart + 1) < size(stimWinFlData(iPos).rawFl, 1)
                currEnd = currEnd + 1;
            end
            stimWinFlData(iPos).rawFl = cat(3, stimWinFlData(iPos).rawFl, ...
                    currTd.rawFlMat(currStart:currEnd, :));     % --> [vol, ROI, stim]
            stimWinFlData(iPos).dff = cat(3, stimWinFlData(iPos).dff, ...
                    currTd.dffMat(currStart:currEnd, :));       % --> [vol, ROI, stim]
            stimWinFlData(iPos).expDff = cat(3, stimWinFlData(iPos).expDff, ...
                    currTd.expDffMat(currStart:currEnd, :));    % --> [vol, ROI, stim]
            stimWinFlData(iPos).zscore = cat(3, stimWinFlData(iPos).zscore, ...
                    currTd.zscoreMat(currStart:currEnd, :));    % --> [vol, ROI, stim]
        end
       
        % Extract optic flow data for all trials from the stim window
        for iStim = 1:numel(stimWinStartFlowFrames)
            stimFlowData = currFlow(stimWinStartFlowFrames(iStim):stimWinEndFlowFrames(iStim));
            if ~isempty(stimWinFlow) && numel(stimFlowData) > size(stimWinFlow, 3)
                stimWinFlow(iTrial, iStim, :) = stimFlowData(1:size(stimWinFlow, 3)); % --> [trial, stim, frame]
            else
                stimWinFlow(iTrial, iStim, :) = stimFlowData;
            end
        end
        
    end%iPos    
    
    
end%iTrial

for iPos = 1:numel(stimWinFlData)
   disp(['Stim position ', num2str(iPos), ': ', num2str(size(stimWinFlData(iPos).rawFl, 3)), ...
            ' presentations without fly movement']); 
end
disp(' ')

catch foldME; rethrow(foldME); end

%% Plot mean flow data within stim window

% Get average optic flow within the stim window
smFlow = repeat_smooth(stimWinFlow, 20, 'dim', 3, 'smwin', 5);
avgFlow = squeeze(mean(mean(smFlow, 2), 1));

% Calculate mean FicTrac data frame rate
meanIFI = [];
for iTrial = 1:numel(mD)
    meanIFI(iTrial) = median(diff(mD(iTrial).ftFrameTimes));
end
ftFrameRate = 1/median(meanIFI);

% Plot mean flow data
figure(2);clf; hold on
plotTimes = ((1:numel(avgFlow)) ./ ftFrameRate) - (1/ftFrameRate) ...
    - quiescenceWinSec(1);
plot(plotTimes, avgFlow, 'linewidth', 2)

% Shade stim period
stimDur = stimOffsetTimes(1) - stimOnsetTimes(1);
plot_stim_shading([0 stimDur])
hold on;

%% Plot visual responses of each ROI (excluding epochs with fly movement)

saveFig = 1;

roiNum = 2;
smWin = 5;
singleTrials = 1;
skipTimes = [-0.1 0.38];
avgWinSize = [1 1];

dataType = 'expDff';

try 
% Adjust figure dimensions for number of subplots
if numel(stimWinFlData) == 1
   figDims = [800, 500];
   mr = 0.05;
   ml = 0.12;
else
   figDims = [numel(stimWinFlData) * 600, 500];
   mr = 0.03;
   ml = 0.06;
end

% Create figure
f = figure(1);clf
f.Color = [1 1 1];
f.Position(3:4) = figDims;

roiName = mD(trialNums(1)).roiData(roiNum).name;
clear ax
for iPos = 1:numel(stimWinFlData)
    plotData = squeeze(stimWinFlData(iPos).(dataType)(:, roiNum, :));% - squeeze(stimWinFlData(iPos).rawFl(:, 3, :));
    plotTimes = ((1:size(plotData, 1)) ./ mD(1).volumeRate) - (1/mD(1).volumeRate) ...
            - quiescenceWinSec(1);

    if ~isempty(skipTimes)
        plotData(plotTimes > skipTimes(1) & plotTimes < skipTimes(2), :) = nan;
        plotTimes(plotTimes > skipTimes(1) & plotTimes < skipTimes(2)) = nan;
    end

    ax(iPos) = subaxis(1, numel(stimWinFlData), iPos, ...
            'mb', 0.14, 'mt', 0.15, 'mr', mr, 'ml', ml);
    hold on
    
    % Plot single trial data
    cm = jet(size(plotData, 2)) * 0.95;
    if singleTrials
        for iTrial = 1:size(plotData, 2)
            plot(plotTimes, smoothdata(plotData(:, iTrial), 'gaussian', smWin, 'omitnan'), ...
                    'color', cm(iTrial, :), 'linewidth', 1)
        end
    end
    
    % Shade +/- SEM
    sd = std(plotData, [], 2, 'omitnan');
    avg = mean(smoothdata(plotData, 'gaussian', smWin, 'omitnan'), 2); 
    sd = sd(~isnan(plotTimes));
    avg = avg(~isnan(plotTimes));
    sem = sd ./ (size(plotData, 2)^0.5);
    jbfill(plotTimes(~isnan(plotTimes)), [avg + sem]', [avg - sem]', 'k', 'k', 1, 0.2);
    
    % Plot mean
    plot(plotTimes, mean(smoothdata(plotData, 'gaussian', smWin, 'omitnan'), 2), ...
            'linewidth', 3, 'color', 'k')
    hold on
    
    % Shade stim period
    stimDur = stimOffsetTimes(1) - stimOnsetTimes(1);
    plot_stim_shading([0 stimDur])
    hold on;
    
    % Add labels to plot
    xlim([plotTimes(1), plotTimes(end)]);
    xlabel('Time from stim onset (s)');
    ax(iPos).FontSize = 14;
    stimNameStr = '';
    for iType = 1:numel(stimPosNames)
        if iType == 1
            stimNameStr = regexprep(stimPosNames{1}, ' ', '-');
        else
            stimNameStr = [stimNameStr, '_', regexprep(stimPosNames{iType}, ' ', '-')];
        end
    end
    if strcmp(dataType, 'expDff')
        yLabStr = 'Full exp. dF/F';
        saveName = ['mean_', roiName, '_', lower(stimNameStr), '_responses_expDff'];
    elseif strcmp(dataType, 'dff')
        yLabStr = 'Trial dF/F';
        saveName = ['mean_', roiName, '_', lower(stimNameStr), '_responses_dff'];
    elseif strcmp(dataType, 'zscore')
        saveName = ['mean_', roiName, '_', lower(stimNameStr), '_responses_zscore'];
        yLabStr = 'z-scored dF/F';
    elseif strcmp(dataType, 'rawFl')
        yLabStr = 'Raw Fl';
        saveName = ['mean_', roiName, '_', lower(stimNameStr), '_responses_rawFl'];
    end
    if singleTrials
        saveName = [saveName, '_singleTrials'];
    end
    if iPos == 1
        ylabel(yLabStr);
    end
    
    titleStr = [stimPosNames{iPos}, ' (n=', num2str(size(plotData, 2)), ')'];
    title(titleStr, 'fontsize', 14)
end

linkaxes(ax);

suptitle([[mD(1).expID, '  -  mean ', roiName, ' visual stim responses'], ...
    ' (no fly movement within plot epoch)']);

% Save figure
if saveFig
   saveDir = fullfile(analysisDir, mD(1).expID);
   save_figure(f, saveDir, saveName);
end


if strcmp([stimPosNames{:}], 'Full field flash')
% Plot mean change in dF/F in pre- and post-stim windows

% avgWinSize = 2;
stimStartVol = find(isnan(plotTimes), 1, 'first');
stimEndVol = find(isnan(plotTimes), 1, 'last');
winSizeVols = round(avgWinSize * mD(1).volumeRate);

preStimVols = (stimStartVol - winSizeVols):(stimStartVol - 1);
postStimVols = (stimEndVol + 1):(stimEndVol + winSizeVols);

preStimMean = mean(plotData(preStimVols, :), 1);
postStimMean = mean(plotData(postStimVols, :), 1);

f = figure(2);clf; 
f.Color = [1 1 1];
f.Position(3:4) = [350 425];
hold on

for iTrial = 1:numel(preStimMean) 
    plot([1 2], [preStimMean(iTrial); postStimMean(iTrial)], '-o', 'linewidth', 1, ...
            'color', cm(iTrial, :));
end
xlim([0.5 2.5])

plot([mean(preStimMean), mean(postStimMean)], '-s', 'linewidth', 3, ...
    'markersize', 10, 'color', 'k', 'markerfacecolor', 'k')

ax = gca;
if strcmp(dataType, 'expDff')
    yLabStr = 'Full exp. dF/F';
    titleStr = 'dF/F';
elseif strcmp(dataType, 'dff')
    titleStr = 'dF/F';
    yLabStr = 'Trial dF/F';
elseif strcmp(dataType, 'zscore')
    titleStr = 'dF/F';
    yLabStr = 'z-scored dF/F';
elseif strcmp(dataType, 'rawFl')
    yLabStr = 'Fl';
else
    yLabStr = '';
end
ylabel(yLabStr)
ax.XTick = [1 2];
ax.XTickLabel = {'Pre-stim', 'Post-stim'};
ax.FontSize = 13;
title({[mD(1).expID, '  -  Mean ', roiName, ' ', lower(titleStr)], [' in ', num2str(avgWinSize), ...
        ' sec win pre and post full field flash'], ...
        '(movement excluded)' }, 'fontsize', 11)
ax.YTickLabel{end} = '';

% Save figure
if strcmp(dataType, 'expDff')
    saveName = ['pre-post-stim_mean_', roiName, '_visual_responses_expDff'];
elseif strcmp(dataType, 'dff')
    saveName = ['pre-post-stim_mean_', roiName, '_visual_responses_dff'];
elseif strcmp(dataType, 'zscore')
    saveName = ['pre-post-stim_mean_', roiName, '_visual_responses_zscore'];
elseif strcmp(dataType, 'rawFl')
    saveName = ['pre-post-stim_mean_', roiName, '_visual_responses_rawFl'];
end
if saveFig
    saveDir = fullfile(analysisDir, mD(1).expID);
    save_figure(f, saveDir, saveName);
end

[~, p] = ttest(preStimMean, postStimMean)
end

catch foldME; rethrow(foldME); end

%% Extract fluorescence data for opto/odor stim epochs only when fly is quiescent

trialNums = [1:3 10:14];

quiescenceWinSec = [2 4];

try
stimWinFlData = struct();
flowData = [];
for iTrial = 1:numel(trialNums)
    
    currTrialInd = find([mD.trialNum] == trialNums(iTrial));
    
    % Identify stim epochs
    stimTiming = mD(currTrialInd).optoStimTiming;
    trialDuration = mD(currTrialInd).trialDuration;
    
    stimOnsetFrames = [];
    stimOffsetFrames = [];
    if iTrial == 1
        stimWinFlData.rawFl = [];
        stimWinFlData.dff = [];
        stimWinFlData.expDff = [];
        stimWinFlData.zscore = [];
    end
    
    stimOnsetTimes = stimTiming(1):sum(stimTiming(2:3)):trialDuration;
    stimOffsetTimes = stimOnsetTimes + stimTiming(2);
    stimOnsetFrames = time2idx(stimOnsetTimes, mD(currTrialInd).ftFrameTimes);
    stimOffsetFrames = time2idx(stimOffsetTimes, mD(currTrialInd).ftFrameTimes);
    
    % Filter out stim epochs that are accompanied by movement
    currFlow = meanFlowMags{currTrialInd};
    currFlow(end) = 0;
    flowFrameDur = median(diff(mD(currTrialInd).ftFrameTimes));
    flowFrameTimes = (1:1:numel(currFlow)) * flowFrameDur; % NOTE: this is only an approximation
        
    stimOnsetFlowFrames = []; stimOffsetFlowFrames = [];
    for iStim = 1:numel(stimOnsetTimes)
        [~, stimOnsetFlowFrames(iStim)] = min(abs(stimOnsetTimes(iStim) - flowFrameTimes));
        [~, stimOffsetFlowFrames(iStim)] = min(abs(stimOffsetTimes(iStim) - flowFrameTimes));
    end
    quiescenceWinFlowFrames = ceil(quiescenceWinSec ./ flowFrameDur);
    stimWinStartFlowFrames = stimOnsetFlowFrames - quiescenceWinFlowFrames(1);
    stimWinEndFlowFrames = stimOffsetFlowFrames + quiescenceWinFlowFrames(2);
    
    stimOnsetTimes(stimWinStartFlowFrames < 1) = [];
    stimOffsetTimes(stimWinStartFlowFrames < 1) = [];
    stimWinEndFlowFrames(stimWinStartFlowFrames < 1) = [];
    stimWinStartFlowFrames(stimWinStartFlowFrames < 1) = [];
    
    stimOnsetTimes(stimWinEndFlowFrames >= ...
            numel(mD(currTrialInd).moveFlowFrameDist)) = [];
    stimOffsetTimes(stimWinEndFlowFrames >= ...
            numel(mD(currTrialInd).moveFlowFrameDist)) = [];
    stimWinStartFlowFrames(stimWinEndFlowFrames >= ...
            numel(mD(currTrialInd).moveFlowFrameDist)) = [];
    stimWinEndFlowFrames(stimWinEndFlowFrames >= ...
            numel(mD(currTrialInd).moveFlowFrameDist)) = [];
    
    % Extract optic flow data within stim window
    minMoveDists = [];
    for iStim = 1:numel(stimOnsetTimes)
        minMoveDists(iStim) = min(mD(currTrialInd).moveFlowFrameDist( ...
            stimWinStartFlowFrames(iStim):stimWinEndFlowFrames(iStim)));
        currWinFrames = stimWinStartFlowFrames(iStim):stimWinEndFlowFrames(iStim);  
        if ~isempty(flowData) && numel(currWinFrames) > size(flowData, 2)
            currWinFrames = currWinFrames(1:size(flowData, 2));
        elseif ~isempty(flowData) && numel(currWinFrames) < size(flowData, 2)
            currWinFrames = [currWinFrames(1) - 1, currWinFrames];
        end
        flowData(end + 1, :) = currFlow(currWinFrames);
    end
    
    stimWinOnsetTimes = stimOnsetTimes(logical(minMoveDists)) - quiescenceWinSec(1);
    stimWinOffsetTimes = stimOffsetTimes(logical(minMoveDists)) + quiescenceWinSec(2);
         
    % Extract fluorescence data from target stim window for each ROI
    stimWinOnsetVols = floor(stimWinOnsetTimes .* mD(currTrialInd).volumeRate);
    stimWinOffsetVols = floor(stimWinOffsetTimes .* mD(currTrialInd).volumeRate);
    
    % Drop stimuli if analysis window extends beyond edge of trial
    stimWinOffsetVols(stimWinOnsetVols < 1) = [];
    stimWinOnsetVols(stimWinOnsetVols < 1) = [];
    stimWinOnsetVols(stimWinOffsetVols > currTd.nVolumes) = [];
    stimWinOffsetVols(stimWinOffsetVols > currTd.nVolumes) = [];
    
    for iStim = 1:numel(stimWinOnsetVols)
        currTd = mD(currTrialInd);
        currStart = stimWinOnsetVols(iStim);
        currEnd = stimWinOffsetVols(iStim);
        if ~isempty(stimWinFlData.rawFl) && ...
                    (currEnd - currStart + 1) > size(stimWinFlData.rawFl, 1)
           currEnd = currEnd - 1;
        elseif ~isempty(stimWinFlData.rawFl) && ...
                    (currEnd - currStart + 1) < size(stimWinFlData.rawFl, 1)
            currEnd = currEnd + 1;
        end
        stimWinFlData.rawFl = cat(3, stimWinFlData.rawFl, ...
                currTd.rawFlMat(currStart:currEnd, :));     % --> [vol, ROI, stim]
        stimWinFlData.dff = cat(3, stimWinFlData.dff, ...
                currTd.dffMat(currStart:currEnd, :));       % --> [vol, ROI, stim]
        stimWinFlData.expDff = cat(3, stimWinFlData.expDff, ...
                currTd.expDffMat(currStart:currEnd, :));    % --> [vol, ROI, stim]
        stimWinFlData.zscore = cat(3, stimWinFlData.zscore, ...
                currTd.zscoreMat(currStart:currEnd, :));    % --> [vol, ROI, stim]
    end
        
end%iTrial


disp([num2str(size(stimWinFlData(iPos).rawFl, 3)), ' presentations without fly movement']);
disp(' ')

catch foldME; rethrow(foldME); end

%% Plot mean flow data within stim window

% Get average optic flow within the stim window
smFlow = repeat_smooth(flowData, 20, 'dim', 2, 'smwin', 6);
avgFlow = squeeze(mean(smFlow, 1));

% Calculate mean FicTrac data frame rate
meanIFI = [];
for iTrial = 1:numel(mD)
    meanIFI(iTrial) = mean(diff(mD(iTrial).ftFrameTimes));
end
ftFrameRate = 1/mean(meanIFI);

% Plot mean flow data
figure(2);clf; hold on
plotTimes = ((1:numel(avgFlow)) ./ ftFrameRate) - (1/ftFrameRate) ...
    - quiescenceWinSec(1);
plot(plotTimes, avgFlow, 'linewidth', 2)

% Shade stim period
stimDur = stimOffsetTimes(1) - stimOnsetTimes(1);
plot_stim_shading([0 stimDur])
hold on;

%% Plot responses of each ROI

saveFig = 1;


roiNum = 2;
smWin = 5;
singleTrials = 0;
skipTimes = [];

dataType = 'expDff';

stimType = 'ACV_neat';

try 

figDims = [800, 500];
mr = 0.05;
ml = 0.12;


% Create figure
f = figure(1);clf
f.Color = [1 1 1];
f.Position(3:4) = figDims;

roiName = mD(trialNums(1)).roiData(roiNum).name;
clear ax

    plotData = squeeze(stimWinFlData.(dataType)(:, roiNum, :));% - squeeze(stimWinFlData(iPos).rawFl(:, 3, :));
    plotTimes = ((1:size(plotData, 1)) ./ mD(1).volumeRate) - (1/mD(1).volumeRate) ...
            - quiescenceWinSec(1);

    if ~isempty(skipTimes)
        plotData(plotTimes > skipTimes(1) & plotTimes < skipTimes(2), :) = nan;
        plotTimes(plotTimes > skipTimes(1) & plotTimes < skipTimes(2)) = nan;
    end

    ax = subaxis(1,1,1, ...
            'mb', 0.14, 'mt', 0.15, 'mr', mr, 'ml', ml);
    hold on
    
    % Plot single trial data
    cm = jet(size(plotData, 2)) * 0.95;
    if singleTrials
        for iTrial = 1:size(plotData, 2)
            plot(plotTimes, smoothdata(plotData(:, iTrial), 'gaussian', smWin, 'omitnan'), ...
                    'color', cm(iTrial, :), 'linewidth', 1)
        end
    end
    
    % Shade +/- SEM
    sd = std(plotData, [], 2, 'omitnan');
    avg = mean(smoothdata(plotData, 'gaussian', smWin, 'omitnan'), 2); 
    sd = sd(~isnan(plotTimes));
    avg = avg(~isnan(plotTimes));
    sem = sd ./ (size(plotData, 2)^0.5);
    jbfill(plotTimes(~isnan(plotTimes)), [avg + sem]', [avg - sem]', 'k', 'k', 1, 0.2);
    
    % Plot mean
    plot(plotTimes, mean(smoothdata(plotData, 'gaussian', smWin, 'omitnan'), 2), ...
            'linewidth', 3, 'color', 'k')
    hold on
    
    % Shade stim period
    stimDur = stimOffsetTimes(1) - stimOnsetTimes(1);
    plot_stim_shading([0 stimDur])
    hold on;
    
    % Add labels to plot
    xlim([plotTimes(1), plotTimes(end)]);
    xlabel('Time from stim onset (s)');
    ax(iPos).FontSize = 14;

    if strcmp(dataType, 'expDff')
        yLabStr = 'Full exp. dF/F';
        saveName = ['mean_', roiName, '_', regexprep(stimType, '_', '\_'), '_responses_expDff'];
    elseif strcmp(dataType, 'dff')
        yLabStr = 'Trial dF/F';
        saveName = ['mean_', roiName, '_', regexprep(stimType, '_', '\_'), '_responses_dff'];
    elseif strcmp(dataType, 'zscore')
        saveName = ['mean_', roiName, '_', regexprep(stimType, '_', '\_'), '_responses_zscore'];
        yLabStr = 'z-scored dF/F';
    elseif strcmp(dataType, 'rawFl')
        yLabStr = 'Raw Fl';
        saveName = ['mean_', roiName, '_', regexprep(stimType, '_', '\_'), '_responses_rawFl'];
    end
    if singleTrials
        saveName = [saveName, '_singleTrials'];
    end
    if iPos == 1
        ylabel(yLabStr);
    end

linkaxes(ax);

title(sprintf('Mean %s responses (n=%d)', regexprep(stimType, '_', '\\_'), size(plotData, 2)))

suptitle([[mD(1).expID, '  - ', roiName, ' responses'], ...
    ' (no fly movement within plot epoch)']);

% Save figure
if saveFig
   saveDir = fullfile(analysisDir, mD(1).expID);
   save_figure(f, saveDir, saveName);
end



catch foldME; rethrow(foldME); end

