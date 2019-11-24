

currTrial = 6;

td = mD([mD.trialNum] == currTrial);

nROIs = numel(td.roiData);


% Create matrix of dF/F data for all ROIs in the current trial
dffMat = []; roiInds = [];
roiNames = {td.roiData.name};
for iGlom = 1:16
    % Rearrange the ROIs if necessary to make sure they're in the right order 
    if iGlom <= 8
        glomName = ['L', num2str(iGlom)];
    else
        glomName = ['R', num2str(16 - iGlom + 1)];
    end
    glomInd = find(ismember(roiNames, glomName));
    if isempty(glomInd)
        dffMat(:, end + 1) = nan(td.nVolumes, 1);
    else
        dffMat(:, end + 1) = td.roiData(glomInd).dffData;
    end
end

% Average across the two sides of the PB
wedgeDffMat = sum(cat(3, dffMat(:, 1:8), dffMat(:, 9:16)), 3, 'omitnan') ./ 2;

% Calculate the population vector average and the amplitude of the summed population vector
% (AKA 'PVA strength') for each volume
meanMat = wedgeDffMat; 


% meanMat(:, 8) = mean([meanMat(:,1), meanMat(:, 7)], 2); % Temp
% meanMat(:, 8) = 0;

angleMat = repmat((-7 * pi/8):pi/4:(7 * pi/8), td.nVolumes, 1);
[x, y] = pol2cart(angleMat, meanMat); % Convert to cartesian coordinates to add vectors
[theta, rho] = cart2pol(sum(x, 2), sum(y, 2));
theta = -theta; % Inverting polarity so it matches other data in plots
dffVectStrength = rho;
dffVectAvgRad = theta; 
dffVectAvgWedge = (theta/pi * 4 + 4.5);



% PLOT DATA THROUGHOUT TRIAL

figure(1);clf;
clear allAx

subplot(511)
allAx(1) = gca;
% dF/F population vector strength
plot(td.volTimes, dffVectStrength, 'color', rgb('darkgreen'), 'linewidth', 1.5);
if td.usingOptoStim
   hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
end
ylabel('PVA strength')
ylim([0 10]);

subplot(512)
allAx(2) = gca;
% Mean dF/F data for each glomerulus
imagesc([0, td.trialDuration], [1, 8], wedgeDffMat');%colorbar
hold on; % Overlay the dF/F population vector average
plot(td.volTimes, max(dffVectAvgWedge) - dffVectAvgWedge, ... % Inverting because of imagesc axes
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
% uwVectAvgRad = unwrap(dffVectAvgRad + pi);
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


%% Calculate average responses at each bar position







