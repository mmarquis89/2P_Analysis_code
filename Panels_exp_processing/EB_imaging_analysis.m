

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
roiData.expBaseline = roiData.baselineVals;

% Load FicTrac data
ftFile = fullfile(dataDir, [expID, '_ficTracData.mat']);
load(ftFile, 'ftData');

% Load flailing event data
eventData = load_event_data({expID}, dataDir);
flailingEvents = eventData.flailing;

% Load panels metadata
panelsMdFile = fullfile(dataDir, [expID, '_panelsMetadata.mat']);
load(panelsMdFile, 'panelsMetadata');

% Load reference images 
refImgFile = fullfile(dataDir, [expID, '_refImages.mat']);
load(refImgFile, 'refImages');


%% Plot summary of movement throughout experiment

allFlow = cell2mat(ftData.meanFlow);

%% Plot ROI data

roiNames = {'BU-L', 'BU-R'};
smWin = 5;

trialNums = 4;

f = figure(2);clf; hold on;
f.Color = [1 1 1];

allFl = cell2mat(roiData(strcmp(roiData.roiName, roiNames{1}), :).rawFl');
ax = subaxis(3, 1, 1);
plot(smoothdata(allFl(:, trialNums), 1, 'gaussian', smWin));

allFl = cell2mat(roiData(strcmp(roiData.roiName, roiNames{2}), :).rawFl');
ax = subaxis(3, 1, 3);
plot((smoothdata(allFl(:, trialNums), 1, 'gaussian', smWin)));

ax = subaxis(3, 1, 2);
% allPosX = cell2mat(panelsMetadata.panelsPosX)
plot(1:numel(panelsMetadata.panelsPosX{1}), as_vector(panelsMetadata.panelsPosX{1}'))

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


