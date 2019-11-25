
currTrial = 6;
td = mD([mD.trialNum] == currTrial);


%% PLOT DATA THROUGHOUT TRIAL

figure(1);clf;
clear allAx

subplot(511)
allAx(1) = gca;
% dF/F population vector strength
plot(td.volTimes, td.dffVectStrength, 'color', rgb('darkgreen'), 'linewidth', 1.5);
if td.usingOptoStim
   hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
end
ylabel('PVA strength')
ylim([0 10]);

subplot(512)
allAx(2) = gca;
% Mean dF/F data for each glomerulus
imagesc([0, td.trialDuration], [1, 8], td.wedgeDffMat');%colorbar
hold on; % Overlay the dF/F population vector average
plot(td.volTimes, max(td.dffVectAvgWedge) - td.dffVectAvgWedge, ... % Inverting because of imagesc axes
        'color', 'r', 'linewidth', 1.25);
ylabel('dF/F and PVA');
    
subplot(513)
allAx(3) = gca;
% Panels bar position
if td.usingPanels
    plot(td.panelsFrameTimes, td.panelsPosX, 'linewidth', 1.5, 'color', 'b');%colorbar
end
if td.usingOptoStim
   hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
end
ylabel('Bar position')
allAx(3).YTickLabel = [];

subplot(514)
allAx(4) = gca;
% FicTrac heading overlaid with dF/F pva
HD = repeat_smooth(unwrap(td.ftData.intHD), 1, 'smWin', 1);
plot(td.ftFrameTimes, mod(HD, 2*pi), 'color', rgb('purple'));
% hold on
% uwVectAvgRad = unwrap(td.dffVectAvgRad + pi);
% uwVectAvgRad = uwVectAvgRad - uwVectAvgRad(1) + 2*pi;
% plot(td.volTimes, mod(uwVectAvgRad, 2*pi));
if td.usingOptoStim
   hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
end
ylim([0 2*pi])
ylabel('Fly heading (rad)')

subplot(515)
allAx(5) = gca;
% FicTrac movement speed
plot(td.ftFrameTimes, repeat_smooth(td.ftData.moveSpeed, 15, 'smWin', 5), 'color', 'k');%colorbar
if td.usingOptoStim
   hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
end
ylabel('Move speed (mm/sec)')
ylim([0 20])
% subplot(616)
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


%% Calculate average responses to each opto stim period
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
        if startVol > 0 && endVol < td.trialDuration
            stimWinDffData(:, :, iStim) = td.wedgeDffMat(startVol:endVol, :); % --> [volume, wedge, stim]
        end
        stimWinDffAvg = mean(stimWinDffData, 3); % --> [volume, wedge];
        
    end
    plotTimes = td.volTimes(1:size(stimWinDffAvg, 1)) - td.volTimes(baselineVols + 1);
    
    figure(1);clf; 
    plot(plotTimes, smoothdata(stimWinDffAvg, 1, 'gaussian', 3)); 
    hold on; 
    plot_stim_shading([0 (td.optoStimOffsetTimes(1) - td.optoStimOnsetTimes(1))])
end


%% Calculate average responses at each bar position
panelsPosVols = [];
for iVol = 1:size(td.wedgeDffMat, 1)
    [~, currVol] = min(abs(td.panelsFrameTimes - td.volTimes(iVol)));
    panelsPosVols(iVol) = td.panelsPosX(currVol);
end
meanDff = [];
for iPos = 1:numel(unique(td.panelsPosX))
    meanDff(iPos, :) = mean(td.wedgeDffMat(panelsPosVols == (iPos - 1), :), 1); % --> [barPos, wedge]    
end

% Shift so that the center position is directly in front of the fly
meanDffShift = [meanDff(92:96, :); meanDff(1:91, :)];
figure(1); clf;
imagesc([-180:3.75:(180 - 3.75)], [1 8], smoothdata(meanDffShift, 1, 'gaussian', 3)');
hold on; plot([0 0], [0 9], 'color', 'r', 'linewidth', 2)
ylim([0.5 8.5])
xlabel('Bar position (degrees from front of fly)');
ylabel('EB wedge');

%% Average fly movement at each bar position
panelsPosVols = [];
for iVol = 1:size(td.wedgeDffMat, 1)
    [~, currVol] = min(abs(td.panelsFrameTimes - td.volTimes(iVol)));
    panelsPosVols(iVol) = td.panelsPosX(currVol);
end






