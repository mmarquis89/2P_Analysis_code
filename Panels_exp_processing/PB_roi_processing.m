%% Set directories and get expList
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData_60D05_7f_bath_DA';
figDir = fullfile(parentDir, 'Figs');

fileList = dir(fullfile(parentDir, '2020*'));
expList = unique(cellfun(@(x) x(1:10), {fileList.name}, 'uniformoutput', 0));

% ----- Load data for the selected experiments -----

% Load metadata 
[expMd, trialMd] = load_metadata(expList, parentDir);

% Load imaging data
roiData = load_roi_data(expList, parentDir);

% Load FicTrac data
ftData = load_ft_data(expList, parentDir);

% Load flailing event data (if any exists)
try
    eventData = load_event_data(expList, parentDir);
    flailingEvents = eventData.flailing;
catch; end

% Load panels metadata
panelsMetadata = load_panels_metadata(expList, parentDir);


%% CHECK CORRELATION BETWEEN GLOMERULI

expList = unique(expMd.expID);

try
for iExp = 1:size(expMd, 1)
    
    currExpID = expMd.expID{iExp};
    currExpRoiData = roiData(strcmp(roiData.expID, currExpID), :);
    
    % Concatenate data from all trials
    flMat = [];
    for iTrial = 1:max(currExpRoiData.trialNum)
        currTrialRoiData = currExpRoiData(currExpRoiData.trialNum == iTrial, :);
        currFl = (cell2mat(currTrialRoiData.rawFl'));
        flMat = [flMat; currFl];
    end
    
    % Calculate correlations between glomeruli
    R = corrcoef(flMat);
    LRCorrMat = R(1:8, 9:16);

    f = figure(iExp); clf;
    f.Color = [1 1 1];
    imagesc(LRCorrMat);
    ax = gca;
    colormap('bluewhitered')
    ax.YTick = 1:8;
    ax.XTick = 1:8;
    ax.YTickLabel = currTrialRoiData.roiName(1:8);
    ax.XTickLabel = currTrialRoiData.roiName(9:16);
    axis square
    title(currExpID)
    
end
catch ME; rethrow(ME); end
    













