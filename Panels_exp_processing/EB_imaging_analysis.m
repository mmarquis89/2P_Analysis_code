

expDir = uigetdir('D:\Dropbox (HMS)\2P Data\Imaging Data\', 'Select an imaging data directory');
dataDir = fullfile(expDir, 'ProcessedData');

% Scan file names in experiment directory to get the expID
imgDataFiles = dir(fullfile(expDir, ['*trial*.tif']));
expID = imgDataFiles(1).name(1:10);

%% Load data for the selected experiment

% Load metadata 
[expMd, trialMd] = load_metadata({expID}, dataDir);

% Load imaging data
roiData = load_roi_data({expID}, dataDir);
% roiData.expBaseline = roiData.baselineVals;

% Load FicTrac data
ftFile = fullfile(dataDir, [expID, '_ficTracData.mat']);
load(ftFile, 'ftData');

% Load flailing event data
try
    eventData = load_event_data({expID}, dataDir);
    flailingEvents = eventData.flailing;
catch; end

% Load panels metadata
panelsMdFile = fullfile(dataDir, [expID, '_panelsMetadata.mat']);
load(panelsMdFile, 'panelsMetadata');

% Load reference images 
refImgFile = fullfile(dataDir, [expID, '_refImages.mat']);
load(refImgFile, 'refImages');


%% Average 2+ ROIs together to create a new ROI

combROIs = {'EB', 'BU-L', 'BU-R'};
newROI = 'EB-DAN';

% combROIs = {'FB-1', 'FB-2'};
% newROI = 'FB-PPM3';

newRoiData = roiData;
for iTrial = 1:max(roiData.trialNum)
    currTrialRoiData = roiData(roiData.trialNum == iTrial, :);
    currTrialCombRoiData = currTrialRoiData(ismember(currTrialRoiData.roiName, combROIs), :);
    newRow = table(currTrialRoiData.expID(1), currTrialRoiData.trialNum(1), {newROI}, ...
            'variableNames', {'expID', 'trialNum', 'roiName'});
    for iRoi = 1:numel(combROIs)
        if iRoi == 1
            newSubRois = currTrialCombRoiData.subROIs{1};
        else
            newSubRois = [newSubRois, currTrialCombRoiData.subROIs{iRoi}];
        end
    end
    newRow.subROIs = newSubRois;
    newRow.rawFl = {mean(cell2mat(currTrialCombRoiData.rawFl'), 2)};
    newRow.trialBaseline = mean(currTrialCombRoiData.trialBaseline);
    newRow.expBaseline = mean(currTrialCombRoiData.expBaseline);
    newRoiData = [newRoiData; newRow];
end



%% 

currExpID = '20201015-2';
roiName = 'EB';

a = MoveSpeedAnalysis({currExpID}, ... 
        ['D:\Dropbox (HMS)\2P Data\Imaging Data\', currExpID, '_75B10_7f\ProcessedData']);
    
a.params.roiName = roiName;
a.params.nHistBinds = 25;
a.params.flType = 'expDff';
a.params.convertSpeedUnits = 0;
a.params.smWinVols = 3;
a.params.smWinFrames = 5;
a.params.nSmoothRepsFrames = 20;
% a.params.frameRate = 1 / median(cellfun(@(x) median(diff(x)), a.sourceDataTable.sourceData.frameTimes));

a.filterDefs.expID = currExpID;
a.filterDefs.roiName = roiName;
a.init_filters();

a = a.analyze();



%% Plot ROI data

% roiNames = {'EB', 'BU-L', 'BU-R'};
% roiNames = {'FB-1', 'FB-2', 'FB-3', 'FB-4', 'FB-5', 'FB-6'};
% roiNames = {'EB-DAN', 'FB-PPM3', 'FB-6'};
roiNames = {'EB-DAN', 'FB'};

smWin = 5;
smReps = 50;

trialNums = 1;

% currRoiData = roiData;
currRoiData = newRoiData;

f = figure(2);clf; hold on;
f.Color = [1 1 1];

ax = subaxis(2, 1, 1, 'sv', 0.04, 'm', 0.05);
hold on;
flMat = [];
for iROI = 1:numel(roiNames)
    currFl = cell2mat(currRoiData(strcmp(currRoiData.roiName, roiNames{iROI}), :).rawFl');
    flMat = [flMat, currFl(:, trialNums)];
end
flMat = smoothdata(flMat, 1, 'gaussian', smWin);
xx = (1:size(flMat, 1)) ./ expMd.volumeRate;
for iROI = 1:numel(roiNames)
    plot(xx, flMat(:, iROI), 'linewidth', 2);
end
legend(roiNames)
xlim([xx(1), xx(end)])

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

% FicTrac
ax = subaxis(2, 1, 2);
moveSpeed = repeat_smooth(ftData.moveSpeed{trialNums}, smReps, 'smWin', smWin);
fwSpeed = repeat_smooth(ftData.fwSpeed{trialNums}, smReps, 'smWin', smWin);
yawSpeed = repeat_smooth(ftData.yawSpeed{trialNums}, smReps, 'smWin', smWin);
sideSpeed = repeat_smooth(ftData.sideSpeed{trialNums}, smReps, 'smWin', smWin);
totalSpeed = (abs(fwSpeed) + abs(yawSpeed) + abs(sideSpeed)) ./ 3;
frameRate = 1/ median(diff(ftData.frameTimes{trialNums}));
xx = (1:numel(moveSpeed))' ./ frameRate;
plot(xx, moveSpeed, 'linewidth', 2); hold on;
xlim([xx(1), xx(end)])
ylim([0 12])
yyaxis right
plot(xx, abs(yawSpeed), 'linewidth', 2);
legend({'Move', 'Yaw'})
% yyaxis right
% plot((1:numel(moveSpeed))/numel(moveSpeed), abs(yawSpeed), 'linewidth', 1);
% plot((1:numel(totalSpeed))/numel(totalSpeed), totalSpeed, 'linewidth', 1);
% plot((1:numel(fwSpeed))/numel(fwSpeed), abs(fwSpeed), 'linewidth', 1);
% linkaxes()

%% Pairwise scatter plots of all ROIs

roiNames = unique(roiData.roiName);

roiNames = roiNames([3, 4, 8]);

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







%%
smWin = 7;
trialNum = 5;


% Calculate volume and panels frame times
volTimes = calc_volTimes(trialMd.nVolumes(1), expMd.volumeRate, trialMd.trialDuration(1), ...
        trialMd.originalTrialCount(1));
panelsFrameTimes = double(1:trialMd.nPanelsFrames(1)) ./ expMd.panelsDisplayRate;

% Get mean data for each ROI
roiList = unique(roiData.roiName);
roiList = roiList(~strcmp(roiList, 'Background'));
allExpDff = [];
for iRoi = 1:numel(roiList)
    
    currData = roiData(strcmp(roiData.roiName, roiList{iRoi}), :);
    currFl = cell2mat(currData.rawFl'); % --> [volume, trial]
    
    % Calculate full exp dF/F
    allExpDff(:, :, iRoi) = currFl ./ currData.expBaseline(1); % --> [volume, trial, ROI]
end

f = figure(3);clf; hold on;
f.Color = [1 1 1];
ax = subaxis(2, 1, 1, 'sv', 0.1);
plot(repmat(volTimes', 1, size(allExpDff, 3)), ...
        smoothdata(squeeze(allExpDff(:, trialNum, :)), 1, 'gaussian', smWin), 'linewidth', 1);
legend(roiList)
ax.FontSize = 14;
xlabel('Time (sec)');
ylabel('dF/F');
title([expID, ' — 19C08 > Brp-short-GCaMP7f — 100 uM KCl wash-in'])

ax = subaxis(2, 1, 2);
plot(1:numel(panelsMetadata.panelsPosX{1}), as_vector(panelsMetadata.panelsPosX{1}'), 'linewidth', ...
        1.5)
ax.FontSize = 14;
ax.YTick = [];
ax.XTick = [];
xlabel('Bar position');

%% Plot bar-position tuning

smWin = 3;
smReps = 5;
roiNums = [3:10];

% Calculate volume and panels frame times
volTimes = calc_volTimes(trialMd.nVolumes(1), expMd.volumeRate, trialMd.trialDuration(1), ...
        trialMd.originalTrialCount(1));
panelsFrameTimes = double(1:trialMd.nPanelsFrames(1)) ./ expMd.panelsDisplayRate;

% Downsample panels position data to volumes
panelsPosVols = [];
for iVol = 1:size(allFl, 1)
    [~, currVol] = min(abs(panelsFrameTimes - volTimes(iVol)));
    panelsPosVols(iVol) = panelsMetadata.panelsPosX{1}(currVol);
end

% Get mean data for each ROI
roiList = unique(roiData.roiName);
roiList = roiList(~strcmp(roiList, 'Background'));
meanDff = [];
for iRoi = 1:numel(roiList)
    
    currData = roiData(strcmp(roiData.roiName, roiList{iRoi}), :);
    currFl = cell2mat(currData.rawFl'); % --> [volume, trial]
    
    % Calculate full exp dF/F
    expDff = currFl ./ currData.expBaseline(1);
    
    for iPos = 1:numel(unique(panelsMetadata.panelsPosX{1}))
        meanDff(iPos, iRoi, :) = ...
            mean(expDff(panelsPosVols == (iPos - 1), :), 1, 'omitnan'); % --> [barPos, ROI, iTrial]
    end    
end


meanDff = smoothdata(meanDff, 1, 'gaussian', 5, 'omitnan');   % --> [barPos, ROI, iTrial]

% Shift data so center of plot is directly in front of the fly
shiftData = cat(1, meanDff(92:96, :, :), meanDff(1:91, :, :));
shiftDataSm = repeat_smooth(shiftData, smReps, 'smWin', smWin);

fullExpMeanDff = mean(shiftDataSm(:, :, 1:end-1), 3);                            % --> [barPos, ROI]

if ~isempty(roiNums)
    fullExpMeanDff = fullExpMeanDff(:, roiNums);
    roiList = roiList(roiNums);
end

% Plot data
plotX = -180:3.75:(180 - 3.75);

f = figure(1);clf;
f.Color = [1 1 1];
plot(repmat(plotX', 1, size(fullExpMeanDff, 2)), fullExpMeanDff, 'linewidth', 1.5);
legend(roiList, 'location', 'nw')
xlabel('Bar position (deg)');
ylabel ('Mean dF/F')
set(gca, 'FontSize', 14);
title([expID, ' — 20A02 > Brp-short-GCaMP7f'])


