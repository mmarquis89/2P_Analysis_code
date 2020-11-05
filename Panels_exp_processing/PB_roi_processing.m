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
    
glomPairNames = table((1:8)', {'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'}', ...
    {'R1', 'R8', 'R7', 'R6', 'R5', 'R4', 'R3', 'R2'}', 'variablenames', ...
    {'wedge', 'leftGlom', 'rightGlom'});

for iExp = 1:size(expMd, 1)
    
    currExpID = expMd.expID{iExp};
    currExpRoiData = roiData(strcmp(roiData.expID, currExpID), :);
    
    % Concatenate data from all trials
    flMat = [];
    for iTrial = 1:max(currExpRoiData.trialNum)
        currTrialRoiData = currExpRoiData(currExpRoiData.trialNum == iTrial, :);
        
        currLeftFl = [];
        currRightFl = [];
        for iWedge = 1:8
            leftData = currTrialRoiData.rawFl(strcmp(currTrialRoiData.roiName, ...
                    glomPairNames.leftGlom{iWedge}));
            rightData = currTrialRoiData.rawFl(strcmp(currTrialRoiData.roiName, ...
                    glomPairNames.rightGlom{iWedge}));
            if ~isempty(leftData)
                currLeftFl(:, iWedge) = leftData{:};
            else
                currLeftFl(:, iWedge) = nan(size(currTrialRoiData.rawFl{1}));
            end
            if ~isempty(rightData)
                currRightFl(:, iWedge) = rightData{:};
            else
                currRightFl(:, iWedge) = nan(size(currTrialRoiData.rawFl{1}));
            end
        end
        currFl = [currLeftFl, currRightFl];
        flMat = [flMat; currFl];
    end
    
    % Calculate correlations between glomeruli
    R = corrcoef(flMat);
    LRCorrMat = R(1:8, 9:16);

    % Plot figure
    f = figure(iExp); clf;
    f.Color = [1 1 1];
    imagesc(LRCorrMat);
    ax = gca;
    colormap('bluewhitered')
    ax.YTick = 1:8;
    ax.XTick = 1:8;
    ax.YTickLabel = glomPairNames.leftGlom;
    ax.XTickLabel = glomPairNames.rightGlom;
    axis square
    title(currExpID)
    
end
catch ME; rethrow(ME); end
    
%% Calculate PVA and the corresponding amplitude for every imaging volume

wedgeData = roiData(~cellfun(@isempty, regexp(roiData.roiName, 'EB-', 'match')), :);


wedgeData.trialDff = cellfun(@(x, y) (x - y) ./ y, wedgeData.rawFl, ...
        mat2cell(wedgeData.trialBaseline, ones(size(wedgeData.trialBaseline))), 'uniformoutput', 0);
wedgeData.expDff = cellfun(@(x, y) (x - y) ./ y, wedgeData.rawFl, ...
        mat2cell(wedgeData.expBaseline, ones(size(wedgeData.expBaseline))), 'uniformoutput', 0);

%% 







