

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

roiNames = {'BU-L', 'EB-6'};

figure(2);clf; hold on;

allFl = cell2mat(roiData(strcmp(roiData.roiName, roiNames{1}), :).rawFl');
ax = subaxis(3, 1, 1);
plot(smoothdata(allFl, 1, 'gaussian', 5));

allFl = cell2mat(roiData(strcmp(roiData.roiName, roiNames{2}), :).rawFl');
ax = subaxis(3, 1, 3);
plot(smoothdata(allFl, 1, 'gaussian', 5));

ax = subaxis(3, 1, 2);
% allPosX = cell2mat(panelsMetadata.panelsPosX)
plot(panelsMetadata.panelsPosX')




