
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

%%


currTrial = 7;

currTrialInd = find([mD.trialNum] == currTrial);

td = mD(currTrialInd);

flowSmWin = 30;
moveThresh = 0.04;
omitMoveVols = 0;
moveBufferDur = 2;
showPVA = 0;

flData = td.wedgeDffMat;
% flData = td.wedgeRawFlMat;
% flData = td.wedgeZscoreMat;
flData = td.wedgeExpDffMat;

% flData = td.dffMat;
% flData = td.rawFlMat;
% flData = td.zscoreMat;
% flData = td.expDffMat;
% 
% flData = flData(:, [1:16]); 


% % 
%     currGoodVols = goodVols(:, currTrialInd);
%     flData(~currGoodVols, :) = nan;


% % Identify epochs of quiescence vs. flailing
% currTrialFlow = mD(currTrialInd).flowData;
% smFlow = repeat_smooth(currTrialFlow, 20, 'dim', 2, 'smwin', flowSmWin);
% smFlow = smFlow - min(smFlow(:));
% moveFrames = smFlow > moveThresh;
% 
% % For each quiescence frame, find the distance to the nearest movement frame
% moveFrameDist = zeros(1, numel(moveFrames));
% moveFrameInds = find(moveFrames);
% for iFrame = 1:numel(moveFrames)
%     moveFrameDist(iFrame) = min(abs(moveFrameInds - iFrame));
% end
% 
% % Convert to volumes
% volFrames = [];
% volTimes = mD(currTrialInd).volTimes;
% frameTimes = mD(currTrialInd).flowFrameTimes;
% for iVol = 1:size(flData, 1) 
%     [~, volFrames(iVol)] = min(abs(volTimes(iVol) - frameTimes));
% end
% moveVolDist = moveFrameDist(volFrames) .* (numel(volTimes) / numel(frameTimes));
% moveBufferDurVols = round(moveBufferDur * td.volumeRate);
% 
% if omitMoveVols
%    flData(moveVolDist < moveBufferDurVols, :) = nan; 
% end


% Get mean panels pos data
panelsPosVols = [];
for iVol = 1:size(flData, 1)
    [~, currVol] = min(abs(td.panelsFrameTimes - td.volTimes(iVol)));
    panelsPosVols(iVol) = td.panelsPosX(currVol);
end
meanFlData = [];
flDataTemp = flData;
for iPos = 1:numel(unique(td.panelsPosX))
    meanFlData(iPos, :) = ...
            mean(flDataTemp(panelsPosVols == (iPos - 1), :), 1, 'omitnan'); % --> [barPos, wedge]    
end
meanFlShift = cat(1, meanFlData(92:96, :), meanFlData(1:91, :));
cm = hsv(size(meanFlShift, 2)) .* 0.9;

% PLOT DATA THROUGHOUT TRIAL
figure(1);clf;
clear allAx


subaxis(5, 3, 1, 'mb', 0.05, 'mt', 0.03, 'ml', 0.08, 'mr', 0.03)
allAx(1) = gca;
% Heatmap of mean visual tuning
plotX = -180:3.75:(180 - 3.75);
imagesc(plotX, [1 size(meanFlShift, 2)], ...
        smoothdata(meanFlShift, 1, 'gaussian', 3)');
hold on; plot([0 0], [0 size(meanFlShift, 2) + 1], 'color', 'k', 'linewidth', 1)
allAx(1).XTick = -180:45:180;
ylim([0.5 size(meanFlShift, 2) + 0.5])
% xlabel('Bar position (deg)', 'fontsize', 10)
ylabel('Visual tuning')
for iWedge = 1:size(cm, 1)
   plot(plotX(1:6), ones(1, 6) .* iWedge, 'color', cm(iWedge, :), 'linewidth', 3)
   plot(plotX(end - 5:end), ones(1, 6) .* iWedge, 'color', cm(iWedge, :), 'linewidth', 3) 
end


subaxis(5, 3, 2)
allAx(2) = gca;
% Plot of mean visual tuning
hold on;
for iWedge = 1:size(meanFlShift, 2)
    plotX = -180:3.75:(180 - 3.75);
    plotX(isnan(meanFlShift(:, iWedge))) = nan;
    plot(plotX, smoothdata(meanFlShift(:, iWedge), 1, 'gaussian', 4, 'omitnan'), ...
            'color', cm(iWedge, :), 'linewidth', 1);
end
allAx(2).XTick = -180:45:180;
plot([0 0], ylim(), 'color', 'k', 'linewidth', 2)
xlim([-180, 180]);


subaxis(5, 3, 3)
allAx(3) = gca;
% Mean visual tuning in polar coordinates
plotX = deg2rad(-180:3.75:(180 - 3.75));
pax = polaraxes(); hold on
pax.Position = allAx(3).Position;
pax.Position = pax.Position - [0.04 0.025 0 0];
pax.Position = pax.Position .* [1 1 1.2 1.2];
allAx(3).Visible = 'off';
for iWedge = 1:size(meanFlShift, 2)
%     plotData = smoothdata(meanFlShift(:, iWedge)- mean(meanFlShift, 2), 1, 'gaussian', 6) ;
    plotData = smoothdata(meanFlShift(:, iWedge), 1, 'gaussian', 6) ;
    plotData = plotData - min(meanFlShift(:));
    polarplot(plotX, plotData, 'color', cm(iWedge, :), 'linewidth', 1);
end
% polarplot(plotX, smoothdata(mean(meanFlShift, 2)-min(meanFlShift(:)), 1, 'gaussian', 6), 'color', 'k', 'linewidth', 3);
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
pax.ThetaTick = [0:45:179, 225:45:359];
pax.ThetaTickLabel = {'0' '+45', '+90' '+135' '-135', '-90', '-45'};
pax.RTick = [];
% pax.FontSize = 12;


subaxis(5, 3, 4:6)
allAx(end + 1) = gca;
% Mean dF/F data for each glomerulus
imagesc([0, td.trialDuration], [1, 8], smoothdata(flData, 1, 'gaussian', 3)');%colorbar
hold on; 
% Overlay the dF/F population vector average
if showPVA
    plot(td.volTimes, 8.5 - td.dffVectAvgWedge, ... % Inverting because of imagesc axes
        'color', 'r', 'linewidth', 1.25);
end
ylabel('dF/F and PVA');
colorbar;
 

subaxis(5, 3, 7:9)
allAx(end + 1) = gca;
% Panels bar position
if td.usingPanels
    plot(td.panelsFrameTimes, td.panelsPosX, 'linewidth', 1.5, 'color', 'b');%colorbar
    barCenteredFrames = find(td.panelsPosX == 44);
    yL = ylim;
    hold on;
    plotX = td.panelsFrameTimes(barCenteredFrames);
    plot([plotX; plotX], repmat(yL', 1, numel(barCenteredFrames)), 'color', ...
            'r');
%     ylim(yL);
end
if td.usingOptoStim
   hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
end
ylabel('Bar position')
allAx(end).YTickLabel = [];
% allAx(end).XTickLabel = [];
cb = colorbar; 
cb.Visible = 'off';


subaxis(5, 3, 10:12)
allAx(end + 1) = gca;
if exist(fullfile(parentDir, 'flowMags.mat'), 'file')
    
        % Optic flow data to identify flailing
        currFlow = mD(currTrialInd).flowData;
        currFlow(end) = 0;
        plotData = repeat_smooth(currFlow, 20, 'dim', 2, 'smwin', flowSmWin);
        plotData2 = repeat_smooth(currFlow, 20, 'dim', 2, 'smwin', 6);
        plotData = plotData - min(plotData); 
        plotData2 = plotData2 - min(plotData2);
        flowFrameDur = median(diff(mD(currTrialInd).ftFrameTimes));
        flowFrameTimes = (1:1:numel(currFlow)) * flowFrameDur; 
        plot(flowFrameTimes, plotData, 'color', 'k', 'linewidth', 1);
        hold on;
        plot(flowFrameTimes, plotData2, 'color', 'g');
        if td.usingOptoStim
           hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
        end
        plot([td.ftFrameTimes(1), td.ftFrameTimes(end)], [moveThresh, moveThresh],...
                'linewidth', 0.5, 'color', 'r');
        ylim([0 0.5])
else
    % FicTrac heading overlaid with dF/F pva
    HD = repeat_smooth(unwrap(td.ftData.intHD), 1, 'smWin', 1);
    plot(td.ftFrameTimes, 2*pi - mod(HD, 2*pi), 'color', 'k');
    if td.usingOptoStim
        hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
    end
    ylim([0 2*pi])
    ylabel('Fly heading (rad)')
    yyaxis('right')
    uwVectAvgRad = smoothdata(unwrap(td.dffVectAvgRad + pi), 1, 'gaussian', 3);
    uwVectAvgRad = uwVectAvgRad - uwVectAvgRad(1) + 2*pi;
    plot(td.volTimes, mod(uwVectAvgRad, 2*pi), 'color', rgb('darkorange'));
    ylabel('PVA')
    ylim([0 2*pi])
end
cb = colorbar; 
cb.Visible = 'off';


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
cb = colorbar; 
cb.Visible = 'off';
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


currFlData = flData;
R = corrcoef(currFlData);
% R = cov(currFlData);

% LRCorrMat = R(1:8, 9:end);
LRCorrMat = R(1:8, :);
LRCorrMat = R;

figure(2); clf; 
imagesc(LRCorrMat); 
axis equal
ax = gca;
colormap('bluewhitered')
ax.YTick = 1:8;
ax.XTick = 1:8;
ax.YTickLabel = {td.roiData(1:8).name};
ax.XTickLabel = {td.roiData(16:-1:9).name};
% 
% ax.YTick = 1:8;
% ax.XTick = 1:16;
% ax.YTickLabel = {td.roiData(1:8).name};
% ax.XTickLabel = {td.roiData.name};
% ax.Color = [0 0 0]
% hold on
% plot(ax, [8.5, 8.5], [8.5 0.5], 'linewidth', 5, 'color', 'k')
% plot(ax, [0.5, 8.5], [0.5 8.5], 'linewidth', 5, 'color', 'g')
% plot(ax, [8.5, 16.5], [0.5 8.5], 'linewidth', 5, 'color', 'g')

% ax.YTick = 1:16;
% ax.XTick = 1:16;
% ax.YTickLabel = {td.roiData.name};
% ax.XTickLabel = {td.roiData.name};
% ax.Color = [0 0 0]
% hold on
% plot(ax, [8.5, 8.5], [16.5 0.5], 'linewidth', 5, 'color', 'k')
% plot(ax, [0.5, 16.5], [8.5, 8.5], 'linewidth', 5, 'color', 'k')
% plot(ax, [0.5, 16.5], [0.5 16.5], 'linewidth', 3, 'color', 'g')

[~, test] = max(LRCorrMat')

