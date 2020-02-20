
startDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';
parentDir = uigetdir(startDir, 'Select an experiment directory');

analysisDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\Analysis';

load(fullfile(parentDir, 'analysis_data.mat'));


load(fullfile(parentDir, 'flowMags.mat'));


%%


currTrial = 8;
td = mD([mD.trialNum] == currTrial);

plotROIs = [];

smWin = 5;

flData = td.dffMat;
% flData = td.rawFlMat;
% flData = td.zscoreMat;
% flData = td.expDffMat;

try
% % Get mean panels pos data
% panelsPosVols = [];
% for iVol = 1:size(flData, 1)
%     [~, currVol] = min(abs(td.panelsFrameTimes - td.volTimes(iVol)));
%     panelsPosVols(iVol) = td.panelsPosX(currVol);
% end
% meanFlData = [];
% for iPos = 1:numel(unique(td.panelsPosX))
%     meanFlData(iPos, :) = ...
%             mean(flData(panelsPosVols == (iPos - 1), :), 1, 'omitnan'); % --> [barPos, wedge]    
% end
% meanFlShift = cat(1, meanFlData(92:96, :), meanFlData(1:91, :));
% cm = hsv(size(meanFlShift, 2)) .* 0.9;

% PLOT DATA THROUGHOUT TRIAL
figure(1);clf;
clear allAx

subaxis(5, 3, 1, 'mb', 0.05, 'mt', 0.03)
allAx(1) = gca;
% % Heatmap of mean visual tuning
% plotX = -180:3.75:(180 - 3.75);
% imagesc(plotX, [1 size(meanFlShift, 2)], ...
%         smoothdata(meanFlShift, 1, 'gaussian', 3)');
% hold on; plot([0 0], [0 size(meanFlShift, 2) + 1], 'color', 'k', 'linewidth', 1)
% allAx(1).XTick = -180:45:180;
% ylim([0.5 size(meanFlShift, 2) + 0.5])
% % xlabel('Bar position (deg)', 'fontsize', 10)
% ylabel('Visual tuning')
% for iWedge = 1:size(cm, 1)
%    plot(plotX(1:6), ones(1, 6) .* iWedge, 'color', cm(iWedge, :), 'linewidth', 3)
%    plot(plotX(end - 5:end), ones(1, 6) .* iWedge, 'color', cm(iWedge, :), 'linewidth', 3) 
% end


subaxis(5, 3, 2)
allAx(2) = gca;
% % Plot of mean visual tuning
% hold on;
% plotX = -180:3.75:(180 - 3.75);
% for iWedge = 1:size(meanFlShift, 2)
%     plot(plotX, smoothdata(meanFlShift(:, iWedge), 1, 'gaussian', 3), 'color', cm(iWedge, :), ...
%             'linewidth', 1);
% end
% allAx(2).XTick = -180:45:180;
% plot([0 0], ylim(), 'color', 'k', 'linewidth', 2)
% xlim([plotX(1), plotX(end)]);


subaxis(5, 3, 3)
allAx(3) = gca;
% % Mean visual tuning in polar coordinates
% plotX = deg2rad(-180:3.75:(180 - 3.75));
% pax = polaraxes(); hold on
% pax.Position = allAx(3).Position;
% pax.Position = pax.Position - [0.04 0.025 0 0];
% pax.Position = pax.Position .* [1 1 1.2 1.2];
% allAx(3).Visible = 'off';
% for iWedge = 1:size(meanFlShift, 2)
%     plotData = smoothdata(meanFlShift(:, iWedge), 1, 'gaussian', 4);
%     polarplot(plotX, plotData, 'color', cm(iWedge, :), 'linewidth', 1);
% end
% pax.ThetaZeroLocation = 'top';
% pax.ThetaDir = 'clockwise';
% pax.ThetaTick = [0:45:179, 225:45:359];
% pax.ThetaTickLabel = {'0' '+45', '+90' '+135' '-135', '-90', '-45'};
% pax.RTick = [];
% pax.FontSize = 12;

subaxis(5, 3, 4:6)
allAx(end + 1) = gca;
% Mean dF/F data for each ROI
if isempty(plotROIs)
   plotROIs = 1:numel(td.roiData);    
end
% plotData = repeat_smooth(flData(:, plotROIs), 10, 'smwin', 3);
plotData = smoothdata(squeeze(flData(:, plotROIs)), 1, 'gaussian', smWin);
plot(repmat(td.volTimes', 1, size(plotData, 2)), plotData, 'linewidth', 1); 
ylim([min(plotData(:)) max(plotData(:))])
xlim([0 td.trialDuration])
roiNames = {td.roiData(plotROIs).name};
legend(roiNames)
ylabel('dF/F');
    
subaxis(5, 3, 7:9)
allAx(end + 1) = gca;
% Panels bar position
if td.usingPanels
    plot(td.panelsFrameTimes, td.panelsPosX, 'linewidth', 1.5, 'color', 'b');%colorbar
    barCenteredFrames = find(td.panelsPosX == 44);
end
if td.usingOptoStim
   hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
end
ylabel('Bar position')
allAx(end).YTickLabel = [];
% allAx(end).XTickLabel = [];

subaxis(5, 3, 10:12)
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
%         plot(plotData, 'color', 'k');
        hold on;
        if td.usingOptoStim
           hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
        end
%         plot([td.ftFrameTimes(1), td.ftFrameTimes(end)], [flowThresh, flowThresh],...
%                 'linewidth', 1, 'color', 'r');
        ylim([0 1.5])
% 
% % FicTrac heading overlaid with dF/F pva
% HD = repeat_smooth(unwrap(td.ftData.intHD), 1, 'smWin', 1);
% plot(td.ftFrameTimes, 2*pi - mod(HD, 2*pi), 'color', 'k');
% if td.usingOptoStim
%    hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
% end
% ylim([0 2*pi])
% ylabel('Fly heading (rad)')
% yyaxis('right')
% uwVectAvgRad = smoothdata(unwrap(td.dffVectAvgRad + pi), 1, 'gaussian', 3);
% uwVectAvgRad = uwVectAvgRad - uwVectAvgRad(1) + 2*pi;
% plot(td.volTimes, mod(uwVectAvgRad, 2*pi), 'color', rgb('darkorange'));
% ylabel('PVA')
% ylim([0 2*pi])
% % plot(td.ftFrameTimes, repeat_smooth(td.ftData.moveSpeed, 15, 'smWin', 7), 'color', rgb('orange'));
% % ylabel('Move speed (mm/sec)')

subaxis(5, 3, 13:15)
allAx(end + 1) = gca;
% FicTrac movement speed
plot(td.ftFrameTimes, repeat_smooth(td.ftData.moveSpeed, 15, 'smWin', 7), 'color', 'k');%colorbar
if td.usingOptoStim
   hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
end
ylabel('Move speed (mm/sec)')
ylim([0 10])
yyaxis('right')
plot(td.ftFrameTimes, abs(repeat_smooth(td.ftData.yawSpeed, 15, 'smWin', 7)), 'color', rgb('orange'));
ylabel('Yaw speed (rad/sec)');
% subaxis(6, 1, 6)
% allAx(6) = gca;
% % FicTrac yaw velocity
% plot(td.ftFrameTimes, repeat_smooth(td.ftData.yawSpeed, 1, 'smWin', 30), 'color', 'k');
% if td.usingOptoStim
%    hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
% end

xlabel('Time (s)');

% Link the X-axis limits across all plots
% linkaxes(allAx(3:end), 'x');
linkaxes(allAx([4 5 6 7]), 'x');
xlim([0 td.trialDuration])


% Calculate average responses to each opto stim period
baselineDur = 2;
postStimDur = 2;

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
            stimWinDffData(:, :, iStim) = td.wedgeDffMat(startVol:endVol, :); % --> [volume, wedge, stim]
        end
        stimWinDffAvg = mean(stimWinDffData, 3); % --> [volume, wedge];
        
    end
    plotTimes = td.volTimes(1:size(stimWinDffAvg, 1)) - td.volTimes(baselineVols + 1);
    
    figure(2);clf;hold on;
    for iWedge = 1:size(stimWinDffAvg, 2)
        plot(plotTimes, smoothdata(stimWinDffAvg(:, iWedge), 1, 'gaussian', 3), ...
                'color', cm(iWedge, :), 'linewidth', 2);
    end
    plot_stim_shading([0 (td.optoStimOffsetTimes(1) - td.optoStimOnsetTimes(1))])
end
catch foldME; rethrow(foldME); end


%% Identify epochs of quiescence vs. flailing

moveThresh = 0.08;

for iTrial = 1:numel(mD)

   currTrialFlow = meanFlowMags{iTrial};
   smFlow = repeat_smooth(currTrialFlow, 20, 'dim', 2, 'smwin', 6);
   smFlow = smFlow - min(smFlow(:));
   moveFrames = smFlow > moveThresh;
   
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

trialNums = [1:4];
panelsStimPositions = 14;%[19 67 43]; % 
stimPosNames = {'Full field flash'}; % {'Left bar', 'Right bar', 'Center bar'}; %
quiescenceWinSec = [8 8];

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
        stimWinEndFlowFrames(stimWinEndFlowFrames > ...
                numel(mD(currTrialInd).moveFlowFrameDist)) = [];
        stimWinStartFlowFrames(stimWinEndFlowFrames > ...
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
        
        
    end%iPos    
    
    
end%iTrial

for iPos = 1:numel(stimWinFlData)
   disp(['Stim position ', num2str(iPos), ': ', num2str(size(stimWinFlData(iPos).rawFl, 3)), ...
            ' presentations without fly movement']); 
end
disp(' ')

catch foldME; rethrow(foldME); end

%%

saveFig = 2;


roiNum = 2;
smWin = 3;
singleTrials = 1;
skipTimes =[-0.05 0.38];
avgWinSize = 1.5;

dataType = 'dff';


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
    if iPos == 1
        ylabel(yLabStr);
    end

%     title({[mD(1).expID, '  -  mean ', roiName, ' ', stimPosNames{iPos}, ' responses'], ...
%         ' (no fly movement within plot epochs)'}, 'fontsize', 14)
    title(stimPosNames{iPos}, 'fontsize', 14)
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

plot([ones(1, numel(preStimMean)); 2 * ones(1, numel(preStimMean))], ...
         [preStimMean; postStimMean], '-o', ...
         'linewidth', 1);
xlim([0.5 2.5])

plot([mean(preStimMean), mean(postStimMean)], '-s', 'linewidth', 3, ...
    'markersize', 10, 'color', 'k', 'markerfacecolor', 'k')

ax = gca;
if strcmp(dataType, 'expDff')
    yLabStr = 'Full exp. dF/F';
elseif strcmp(dataType, 'dff')
    yLabStr = 'Trial dF/F';
elseif strcmp(dataType, 'zscore')
    yLabStr = 'z-scored dF/F';
elseif strcmp(dataType, 'rawFl')
    yLabStr = 'Raw Fl';
else
    yLabStr = '';
end
ylabel(yLabStr)
ax.XTick = [1 2];
ax.XTickLabel = {'Pre-stim', 'Post-stim'};
ax.FontSize = 14;
title({[mD(1).expID, '  -  Mean ', roiName, ' ', lower(yLabStr)], [' in ', num2str(avgWinSize), ...
        ' sec win pre and post full field flash'], ...
        '(movement excluded)' }, 'fontsize', 11)
ax.YTickLabel{end} = '';

% Save figure
if strcmp(dataType, 'expDff')
    saveName = ['pre-post-stim_mean_', roiName, '_full-field-flash_responses_expDff'];
elseif strcmp(dataType, 'dff')
    saveName = ['pre-post-stim_mean_', roiName, '_full-field-flash_responses_dff'];
elseif strcmp(dataType, 'zscore')
    saveName = ['pre-post-stim_mean_', roiName, '_full-field-flash_responses_zscore'];
elseif strcmp(dataType, 'rawFl')
    saveName = ['pre-post-stim_mean_', roiName, '_full-field-flash_responses_rawFl'];
end
if saveFig
    saveDir = fullfile(analysisDir, mD(1).expID);
    save_figure(f, saveDir, saveName);
end

[~, p] = ttest(preStimMean, postStimMean)
end