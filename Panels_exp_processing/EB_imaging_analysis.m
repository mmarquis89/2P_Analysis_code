%% Set expList

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\EB-DAN_GroupedAnalysisData';
expList = {'20201023-1'};

fileList = dir(fullfile(parentDir, '2020*'));
expList = unique(cellfun(@(x) x(1:10), {fileList.name}, 'uniformoutput', 0));

%% Load data for the selected experiments

% Load metadata 
[expMd, trialMd] = load_metadata(expList, parentDir);

% Load imaging data
roiData = load_roi_data(expList, parentDir);

% Load FicTrac data
ftData = load_ft_data(expList, parentDir);

% Load flailing event data
try
    eventData = load_event_data(expList, parentDir);
    flailingEvents = eventData.flailing;
catch; end

% Load panels metadata
panelsMetadata = load_panels_metadata(expList, parentDir);

%% 

currExpID = '20201023-1';
roiName = 'EB-DAN';
skipTrials = [7:10];

a = MoveSpeedAnalysis({currExpID}, ... 
        ['D:\Dropbox (HMS)\2P Data\Imaging Data\', currExpID, '_75B10_7f\ProcessedData']);
    
a.params.roiName = roiName;
a.params.nHistBinds = 25;
a.params.flType = 'expDff';
a.params.convertSpeedUnits = 0;
a.params.smWinVols = 5;
a.params.smWinFrames = 7;
a.params.nSmoothRepsFrames = 20;
a.params.skipTrials = skipTrials;
% a.params.frameRate = 1 / median(cellfun(@(x) median(diff(x)), a.sourceDataTable.sourceData.frameTimes));

a.filterDefs.expID = currExpID;
a.filterDefs.roiName = roiName;
a = a.init_filters();

a = a.analyze();

%%

speedType = 'move';

a = a.generate_binned_flData([], speedType);
a.plot_binned_fl(gca);




%% Plot ROI data


expID = '20201023-1';

roiNames = {'EB-DAN'};
roiNames = {'EB-DAN'};

smWin = 7;
smReps = 50;

trialNums = 8;

f = figure(2);clf; hold on;
f.Color = [1 1 1];

% Extract ROI data
currRoiData = roiData(strcmp(roiData.expID, expID), :);
flMat = [];
for iROI = 1:numel(roiNames)
    currFl = cell2mat(currRoiData(strcmp(currRoiData.roiName, roiNames{iROI}), :).rawFl');
    flMat = [flMat, currFl(:, trialNums)];
end
flMat = smoothdata(flMat, 1, 'gaussian', smWin);

% Extract FicTrac data
currFtData = ftData(strcmp(ftData.expID, expID), :);
moveSpeed = repeat_smooth(currFtData.moveSpeed{trialNums}, smReps, 'smWin', smWin);
fwSpeed = repeat_smooth(currFtData.fwSpeed{trialNums}, smReps, 'smWin', smWin);
yawSpeed = repeat_smooth(currFtData.yawSpeed{trialNums}, smReps, 'smWin', smWin);
sideSpeed = repeat_smooth(currFtData.sideSpeed{trialNums}, smReps, 'smWin', smWin);
totalSpeed = (abs(fwSpeed) + abs(yawSpeed) + abs(sideSpeed)) ./ 3;
frameRate = 1/ median(diff(currFtData.frameTimes{trialNums}));

% Get downsampled FicTrac data
volTimes = (1:size(flMat, 1)) ./ expMd.volumeRate;
moveSpeedVols = []; fwSpeedVols = []; yawSpeedVols = [];
for iVol = 1:numel(volTimes)
    dsFrame = argmin(abs(currFtData.frameTimes{trialNums} - volTimes(iVol)));
    moveSpeedVols(iVol) = moveSpeed(dsFrame);
    fwSpeedVols(iVol) = fwSpeed(dsFrame);
    yawSpeedVols(iVol) = yawSpeed(dsFrame);
end

flMat = flMat - min(flMat);

% Plot 1
ax1 = subaxis(2, 1, 1, 'sv', 0.04, 'm', 0.05);
hold on;
xx = (1:size(flMat, 1)) ./ expMd(strcmp(expMd.expID, expID), :).volumeRate;
for iROI = 1:numel(roiNames)
    plot(xx, flMat(:, iROI) ./ max(flMat(:, iROI)), 'linewidth', 2, 'color', 'k');
end

xlim([xx(1), xx(end)])

% yyaxis right
plot(xx, fwSpeedVols ./ max(fwSpeedVols), 'linewidth', 2, 'color', 'b')
hold on;
plotYaw = smoothdata(abs(yawSpeedVols), 2, 'gaussian', smWin);
plot(xx, plotYaw ./ max(plotYaw), '-', 'linewidth', 2, 'color', 'r')
legend([roiNames, {'fwSpeed', 'yawSpeed'}])

% % ROI 1 + ROI 2
% allFl = cell2mat(roiData(strcmp(roiData.roiName, roiNames{1}), :).rawFl');
% ax = subaxis(3, 1, 1, 'sv', 0.04, 'm', 0.05);
% plot((1:numel(allFl(:, trialNums)))/numel(allFl(:, trialNums)), smoothdata(allFl(:, trialNums), 1, 'gaussian', smWin), 'linewidth', 1);
% hold on;
% % yyaxis right
% allFl = cell2mat(roiData(strcmp(roiData.roiName, roiNames{2}), :).rawFl');
% plot((1:numel(allFl(:, trialNums)))/numel(allFl(:, trialNums)), smoothdata(allFl(:, trialNums), 1, 'gaussian', smWin), 'linewidth', 1);
% 
% % ROI 3
% allFl = cell2mat(roiData(strcmp(roiData.roiName, roiNames{3}), :).rawFl');
% ax = subaxis(3, 1, 3);
% plot(smoothdata(allFl(:, trialNums), 1, 'gaussian', smWin), 'linewidth', 1);
% 
% % % Panels
% % ax = subaxis(3, 1, 2);
% % % allPosX = cell2mat(panelsMetadata.panelsPosX)
% % plot(1:numel(panelsMetadata.panelsPosX{1}), as_vector(panelsMetadata.panelsPosX{1}'), ...
% %         'linewidth', 1)
% %   

% Plot 2
ax2 = subaxis(2, 1, 2);
xx = (1:numel(moveSpeed))' ./ frameRate;
plot(xx, fwSpeed, 'linewidth', 2); hold on;
xlim([xx(1), xx(end)])
ylim([0 20])
yyaxis right
plot(xx, abs(rad2deg(yawSpeed)), 'linewidth', 2);
legend({'Fw', 'Yaw'})
% yyaxis right
% plot((1:numel(moveSpeed))/numel(moveSpeed), abs(yawSpeed), 'linewidth', 1);
% plot((1:numel(totalSpeed))/numel(totalSpeed), totalSpeed, 'linewidth', 1);
% plot((1:numel(fwSpeed))/numel(fwSpeed), abs(fwSpeed), 'linewidth', 1);
linkaxes([ax1, ax2], 'x')
%

% Calculate fl-speed correlations for each ROI
for iRoi = 1:size(flMat, 2)
    moveR = corrcoef(flMat(:, iRoi), moveSpeedVols');
    fwR = corrcoef(flMat(:, iRoi), abs(fwSpeedVols)');
    yawR = corrcoef(flMat(:, iRoi), abs(yawSpeedVols)');
    
    disp('  ')
    disp(roiNames{iRoi})
    disp(['move: ', num2str(moveR(2,1), 2)])
    disp(['fw: ', num2str(fwR(2,1), 2)])
    disp(['yaw: ', num2str(yawR(2,1), 2)])
end

%% Pairwise scatter plots of all ROIs

roiNames = unique(roiData.roiName);

roiNames = roiNames([4:12]);

nROIs = numel(roiNames);
f = figure(1);clf; hold on;
f.Color = [1 1 1];
ax = subaxis(nROIs, nROIs, 1, 'sv', 0.04, 'm', 0.05);
count = 1;
for i = 1:nROIs
    ROI_1 = roiNames{i};
    fl_1 = cell2mat(roiData(strcmp(roiData.roiName, ROI_1), :).rawFl);
    for j = 1:nROIs
        ROI_2 = roiNames{j};
        fl_2 = cell2mat(roiData(strcmp(roiData.roiName, ROI_2), :).rawFl);
        if i ~= j 
            subaxis(nROIs, nROIs, count)
            scatter(fl_1, fl_2);
            title([ROI_1, ' vs ' ROI_2])
        end
        count = count + 1;
    end
end

