
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\20191119-1_38A11_ChR_60D05_7f\ProcessedData';
analysisDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\Analysis';

load(fullfile(parentDir, 'analysis_data.mat'));

%%
currTrial = 1;

td = mD([mD.trialNum] == currTrial);

flData = td.wedgeDffMat;
% flData = td.wedgeRawFlMat;
% flData = td.wedgeZscoreMat;

% PLOT DATA THROUGHOUT TRIAL
figure(1);clf;
clear allAx

subaxis(5, 1, 1)
allAx(1) = gca;
% dF/F population vector strength
plot(td.volTimes, td.dffVectStrength, 'color', rgb('darkgreen'), 'linewidth', 1.5);
if td.usingOptoStim
   hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
end
ylabel('PVA strength')
ylim([0 10]);

subaxis(5, 1, 2)
allAx(2) = gca;
% Mean dF/F data for each glomerulus
imagesc([0, td.trialDuration], [1, 8], flData');%colorbar
hold on; % Overlay the dF/F population vector average
plot(td.volTimes, max(td.dffVectAvgWedge) - td.dffVectAvgWedge, ... % Inverting because of imagesc axes
        'color', 'r', 'linewidth', 1.25);
ylabel('dF/F and PVA');
    
subaxis(5, 1, 3)
allAx(3) = gca;
% Panels bar position
if td.usingPanels
    plot(td.panelsFrameTimes, td.panelsPosX, 'linewidth', 1.5, 'color', 'b');%colorbar
    barCenteredFrames = find(td.panelsPosX == 44);
    yL = ylim;
    hold on;
    plotX = td.panelsFrameTimes(barCenteredFrames);
    plot([plotX; plotX], repmat(yL', 1, numel(barCenteredFrames)), 'color', ...
            'r');
%T     ylim(yL);
end
if td.usingOptoStim
   hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
end
ylabel('Bar position')
allAx(3).YTickLabel = [];

subaxis(5, 1, 4)
allAx(4) = gca;
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

subaxis(5, 1, 5)
allAx(5) = gca;
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
linkaxes(allAx, 'x');
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
