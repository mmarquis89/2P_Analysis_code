parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\20200203-4_38A11_ChR_60D05_7f\ProcessedData';

% Load left PB ROI data
load(fullfile(parentDir, 'roiData_reg_L.mat'), 'allROIData');
leftROIData = allROIData;
leftDffMat = cell2mat({leftROIData.roiDefs.dffData}');
leftROINames = {leftROIData.roiDefs.name};

% Load right PB ROI data
load(fullfile(parentDir, 'roiData_reg_R.mat'), 'allROIData');
rightROIData = allROIData;
rightDffMat = cell2mat({rightROIData.roiDefs.dffData}');
rightROINames = {rightROIData.roiDefs.name};

fullDffMat = [leftDffMat', rightDffMat'];
allROINames = [leftROINames, rightROINames];

allR = corrcoef(fullDffMat);

leftR = allR(1:numel(leftROINames), 1:numel(leftROINames));
rightR = allR(numel(leftROINames) + 1:end, numel(leftROINames) + 1:end);

%%

figure(1); clf
imagesc(allR);
axis square
ax = gca;
colormap('bluewhitered')

%%
figure(1);clf
% imagesc(allR(1:15, [1:15, 32 31 29 30 28:-1:16]))
imagesc(allR(end:-1:16, [1:17, 32:-1:18]))
hold on
plot([17.5 17.5],[0.5 17.5],   'color', 'k', 'linewidth', 3)
axis equal
ax = gca;
ax.YLim = [0.5 17.5];

%% Clustering

clusterData = fullDffMat;

[idx, C] = kmeans(clusterData, 8);

[~, idx_test] = pdist2(C, rightDffMat, 'euclidean', 'Smallest', 1);

test = clusterdata(clusterData, 8)



%%
currFlData = leftDffMat;
R = corrcoef(currFlData);
% R = cov(currFlData);
LRCorrMat = R(1:8, 9:end);
% LRCorrMat = R;
figure(2); clf; 
imagesc(LRCorrMat); 
axis square
ax = gca;
colormap('bluewhitered')

ax.YTick = 1:8;
ax.XTick = 1:8;
ax.YTickLabel = {td.roiData(1:8).name};
ax.XTickLabel = {td.roiData(16:-1:9).name};

% ax.YTick = 1:16;
% ax.XTick = 1:16;
% ax.YTickLabel = {td.roiData([1:8, 16:-1:9]).name};
% ax.XTickLabel = {td.roiData([1:8, 16:-1:9]).name};

[~, test] = max(LRCorrMat')