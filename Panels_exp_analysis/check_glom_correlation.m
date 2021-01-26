function LRCorrMat = check_glom_correlation(expID, expMd, roiData)
% ==================================================================================================
% Calculates and plots the correlation coefficients between 


% ==================================================================================================

glomPairNames = table((1:8)', {'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'}', ...
    {'R1', 'R8', 'R7', 'R6', 'R5', 'R4', 'R3', 'R2'}', 'variablenames', ...
    {'wedge', 'leftGlom', 'rightGlom'});

% Find index of current expID (for choosing fig number)
fullExpList = unique(expMd.expID);
currExpInd = find(strcmp(fullExpList, expID), 1);

% Return empty if expID is not in expMD
if isempty(currExpInd)
    LRCorrMat = [];
    return
end

currExpRoiData = roiData(strcmp(roiData.expID, expID), :);

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
end%iTrial

% Calculate correlations between glomeruli
R = corrcoef(flMat);
LRCorrMat = R(1:8, 9:16);

% Plot figure
f = figure(currExpInd); clf;
f.Color = [1 1 1];
imagesc(LRCorrMat);
ax = gca;
colormap('bluewhitered')
ax.YTick = 1:8;
ax.XTick = 1:8;
ax.YTickLabel = glomPairNames.leftGlom;
ax.XTickLabel = glomPairNames.rightGlom;
axis square
title(expID)

end