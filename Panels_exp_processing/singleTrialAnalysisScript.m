
startDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';
parentDir = uigetdir(startDir, 'Select an experiment directory');

analysisDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\Analysis';

load(fullfile(parentDir, 'analysis_data.mat'));


% load(fullfile(parentDir, 'flowMags.mat'));


%%


currTrial = 5

td = mD([mD.trialNum] == currTrial);

flData = td.wedgeDffMat;
% flData = td.wedgeRawFlMat;
% flData = td.wedgeZscoreMat;

% flData = td.dffMat;
% flData = td.rawFlMat;


% Get mean panels pos data
panelsPosVols = [];
for iVol = 1:size(flData, 1)
    [~, currVol] = min(abs(td.panelsFrameTimes - td.volTimes(iVol)));
    panelsPosVols(iVol) = td.panelsPosX(currVol);
end
meanFlData = [];
for iPos = 1:numel(unique(td.panelsPosX))
    meanFlData(iPos, :) = ...
            mean(flData(panelsPosVols == (iPos - 1), :), 1); % --> [barPos, wedge]    
end
meanFlShift = cat(1, meanFlData(92:96, :), meanFlData(1:91, :));
cm = hsv(size(meanFlShift, 2)) .* 0.9;

% PLOT DATA THROUGHOUT TRIAL
figure(1);clf;
clear allAx

subaxis(5, 3, 1, 'mb', 0.05, 'mt', 0.03)
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
plotX = -180:3.75:(180 - 3.75);
for iWedge = 1:size(meanFlShift, 2)
    plot(plotX, smoothdata(meanFlShift(:, iWedge), 1, 'gaussian', 3), 'color', cm(iWedge, :), ...
            'linewidth', 1);
end
allAx(2).XTick = -180:45:180;
plot([0 0], ylim(), 'color', 'k', 'linewidth', 2)
xlim([plotX(1), plotX(end)]);


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
    plotData = smoothdata(meanFlShift(:, iWedge), 1, 'gaussian', 4);
    polarplot(plotX, plotData, 'color', cm(iWedge, :), 'linewidth', 1);
end
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
hold on; % Overlay the dF/F population vector average
plot(td.volTimes, max(td.dffVectAvgWedge) - td.dffVectAvgWedge, ... % Inverting because of imagesc axes
        'color', 'r', 'linewidth', 1.25);
ylabel('dF/F and PVA');
    
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

subaxis(5, 3, 10:12)
allAx(end + 1) = gca;
%         % Optic flow data to identify flailing
%         flowThresh = 0.08;
%         currFlow = meanFlowMags{[mD.trialNum] == currTrial};
%         currFlow(end) = 0;
%         plotData = repeat_smooth(currFlow, 20, 'dim', 2, 'smwin', 6);
%         plotData = plotData - min(plotData); 
% %         plotData(plotData < flowThresh) = nan;
% %         plotData(end) = nan;
%         plot(td.ftFrameTimes(1:numel(plotData)), plotData, 'color', 'k');
%         hold on;
%         if td.usingOptoStim
%            hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
%         end
%         plot([td.ftFrameTimes(1), td.ftFrameTimes(end)], [flowThresh, flowThresh],...
%                 'linewidth', 1, 'color', 'r');
%         ylim([0 0.5])


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
% plot(td.ftFrameTimes, repeat_smooth(td.ftData.moveSpeed, 15, 'smWin', 7), 'color', rgb('orange'));
% ylabel('Move speed (mm/sec)')

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
linkaxes(allAx(3:end), 'x');
xlim([0 td.trialDuration])


% Calculate average responses to each opto stim period
baselineDur = 1;
postStimDur = 1;

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
    
    figure(2);clf; 
    plot(plotTimes, smoothdata(stimWinDffAvg, 1, 'gaussian', 3)); 
    hold on; 
    plot_stim_shading([0 (td.optoStimOffsetTimes(1) - td.optoStimOnsetTimes(1))])
end
