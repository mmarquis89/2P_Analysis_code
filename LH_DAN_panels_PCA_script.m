

startDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';
parentDir = uigetdir(startDir, 'Select an experiment directory');


%% Run PCA on entire dataset

trialNum = [];

imgDataFiles = dir(fullfile(parentDir, 'imagingData_reg_trial_*.mat'));

nPCs = 30;

try 
    
if ~isempty(trialNum)
    imgDataFiles = imgDataFiles(trialNum);
end

for iFile = 1:numel(imgDataFiles)
    
    % Load registered data for current trial
    load(fullfile(parentDir, ['imagingData_reg_trial_', ...
            pad(num2str(iFile), 3, 'left', '0'), '.mat']), 'imgData');
    disp(['Trial #', num2str(iFile), ' loaded'])
    
    sz = size(imgData);
    
    % Pre-smooth with gaussian filter
    for iVol = 1:sz(4)
        for iPlane = 1:sz(3)
            imgData(:, :, iPlane, iVol) = imgaussfilt(imgData(:, :, iPlane, iVol), 1);
        end
    end
    disp('Smoothing complete')
    
    % Run PCA
    data2D = double(reshape(imgData, [(sz(1)*sz(2)*sz(3)), sz(4)]));     % --> [pixel, trialVolume]
    clear imgData
    tData2D = data2D';                                                   % --> [trialVolume, pixel]
    clear data2D
    tic
    [coeff,score,latent,~,explained] = pca(tData2D, 'numcomponents', nPCs); % [pixel, pc]
    tc = toc;
    pcaData = reshape(coeff, [sz(1:3), nPCs]);                        % [y, x, plane, pc]
    clear tData2D
    
    fileName = ['pcaData_trial_', pad(num2str(iFile), 3, 'left', '0'), '.mat'];
    save(fullfile(parentDir, fileName), 'pcaData')
    disp(['Trial ', num2str(iFile), ' PCA completed in ', num2str(tc), ' sec'])
end

catch ME; rethrow(ME); end

%% Run PCA restricted to specific ROIs within each plane

trialNums = [1 2 6 7];

imgDataFiles = dir(fullfile(parentDir, 'imagingData_reg_trial_*.mat'));

nPCs = 30;

try 
    
if ~isempty(trialNums)
    imgDataFiles = imgDataFiles(trialNums);
end

for iFile = 1:numel(imgDataFiles)
    
    % Load registered data for current trial
    disp('Loading imaging data...')
    load(fullfile(parentDir, imgDataFiles(iFile).name), 'imgData');
    trialNumStr = regexp(imgDataFiles(iFile).name, '0..(?=.mat)', 'match');
    trialNumStr = trialNumStr{:};
    
    sz = size(imgData);
    
    % Smooth raw input data
    disp('Smoothing imaging data...')
    for iVol = 1:sz(4)
        for iPlane = 1:sz(3)
            imgData(:, :, iPlane, iVol) = imgaussfilt(imgData(:, :, iPlane, iVol), 1);
        end
    end
    
    % Load ROI
    disp('Loading and applying ROI data...')
    load(fullfile(parentDir, ['roiDefs_PCA_trial_', trialNumStr, '.mat']), 'roiDefs');
    planeList = unique([roiDefs.subROIs.plane]);
    for iPlane = 1:numel(planeList)
        currPlane = planeList(iPlane);
        currSubRois = roiDefs.subROIs([roiDefs.subROIs.plane] == currPlane);
        currPlaneData = squeeze(imgData(:, :, currPlane, :));
        currPlaneMask = any(reshape([currSubRois.mask], size(imgData, 1), size(imgData, 2), []), 3);
        currPlaneMask = repmat(currPlaneMask, 1, 1, size(currPlaneData, 3));
        currPlaneData(~currPlaneMask) = nan;
        imgData(:, :, currPlane, :) = currPlaneData;
    end
    
    % Run PCA on data within ROIs
    disp('Running PCA...')
    data2D = double(reshape(imgData, [(sz(1)*sz(2)*sz(3)), sz(4)])');     % --> [trialVolume, pixel]
    clear imgData
    dataCols = ~isnan(data2D(1, :));
    tic
    [coeff,score,latent,~,explained] = pca(data2D(:, dataCols), 'numcomponents', nPCs); % [pixel, pc]
    tc = toc;
    
    % Return output coeffs to original shape
    pcaData = nan(size(data2D, 2), nPCs);
    pcaData(repmat(dataCols', 1, nPCs)) = coeff;
    pcaData = reshape(pcaData, [sz(1:3), nPCs]);                        % [y, x, plane, pc]
    
    fileName = ['pcaData_ROI_trial_', trialNumStr, '.mat'];
    save(fullfile(parentDir, fileName), 'pcaData')
    disp(['Trial ', trialNumStr, ' PCA completed in ', num2str(tc), ' sec'])
end


catch ME; rethrow(ME); end

%% Load and plot PCA data

% currTrials = [1 2 3 11 12 13]

trialNum = 6; 
% trialNum = currTrials(6)

nPCs = 9;
offset = 0;
targetPCs = [3 7];
useTargetPCs = 0;
sigma = 0.4;
maxIntensity = 400;
pcCLimScaleFactor = 0.4;

load(fullfile(parentDir, ['pcaData_ROI_trial_', pad(num2str(trialNum), 3, 'left', '0'), '.mat']), ...
        'pcaData')
if sum(isnan(pcaData(:))) == 0
   pcaData(pcaData == 0) = nan; 
end
    
% load(fullfile(parentDir, ['pcaData_trial_', pad(num2str(trialNum), 3, 'left', '0'), '.mat']), ...
%         'pcaData')

try 
load(fullfile(parentDir, ['refImages_reg_trial_', ...
        pad(num2str(trialNum), 3, 'left', '0'), '.mat']), 'refImages')
    
    
if ~exist('trialPCs', 'var')
    trialPCs = {};
end
if ~isempty(targetPCs)
    trialPCs{trialNum} = targetPCs;
elseif useTargetPCs
    targetPCs = trialPCs{trialNum}
end

threshVals = [0 0 0];
% threshVals = [0.001, 0.001, .002];

if ~isempty(targetPCs)
    nPCs = numel(targetPCs);
end
    
nPlanes = size(pcaData, 3);
f = figure(2); clf;
f.Color = [1 1 1];
subaxis(nPCs + 1, nPlanes, 1, 1, ...
        'm', 0.01, ...
        'p', 0, ...
        'sh', 0.001, 'sv', 0)

for iPlane = 1:nPlanes
    
    % Plot reference images
    subaxis(nPCs + 1, nPlanes, iPlane, 1)
    imshow(refImages(:, :, iPlane), [0 maxIntensity]);
    colormap(gca, 'gray');
    
    % Plot top PCs
    for iPC = 1:nPCs
       
       if isempty(targetPCs)
           currData = squeeze(pcaData(:, :, iPlane, iPC + offset));
       else
           currData = squeeze(pcaData(:, :, iPlane, targetPCs(iPC)));
       end
       currData = flipud(currData);
       
       % Plot reference image underneath
       ax1 = subaxis(nPCs + 1, nPlanes, iPlane, iPC + 1); hold on
       h1 = imagesc(ax1, flipud(refImages(:, :, iPlane))); axis image
       ax1.CLim = [0 maxIntensity];
       colormap(ax1, 'gray')
       
       ax2 = axes; hold on;
       ax2.Position = ax1.Position;
       
       currData = imgaussfilt(currData, sigma);
       
       h2 = imagesc(ax2, currData, 'alphadata', ~isnan(currData));
       colormap(ax2, 'bluewhitered'); axis image
       ax2.CLim = ax2.CLim * pcCLimScaleFactor;
%        colorbar
       
       ax1.Visible = 'off';
       ax2.Visible = 'off';
       
    end
end

catch ME; rethrow(ME); end 

%% Compare two PCs by overlaying thresholded versions

planeNum = 5;

maxIntensity = 300;
compPCs = [3 7];
directions = [1  -1]; 
thresholds = [30 30 12];
colors = {'magenta', 'lime', 'yellow'};

pc_1 = 10000 * squeeze(pcaData(:, :, planeNum, compPCs(1))) .* directions(1);
pc_2 = 10000 * squeeze(pcaData(:, :, planeNum, compPCs(2))) .* directions(2);

pc_1_thresh = pc_1 > thresholds(1);
pc_2_thresh = pc_2 > thresholds(2);

if numel(compPCs) > 2
    pc_3 = 10000 * squeeze(pcaData(:, :, planeNum, compPCs(3))) .* directions(3);
    pc_3_thresh = pc_3 > thresholds(3);
end

combData = pc_1_thresh + (2 * pc_2_thresh);
if numel(compPCs) > 2
    combData = combData + (4 * pc_3_thresh);
end
combData(combData == 0) = nan;

f = figure(7);clf; 
f.Color = [1 1 1];

% Plot reference image in top subplot
subaxis(2, 1, 1, 'ml', 0.03, 'mr', 0.03, 'mb', 0.03, 'mt', 0.03)
imshow(refImages(:, :, planeNum), [0 maxIntensity]);
colormap(gca, 'gray');
% 
% Plot reference image underneath bottom subplot
ax1 = subaxis(2, 1, 2); hold on
h1 = imagesc(ax1, flipud(refImages(:, :, planeNum))); axis image
ax1.CLim = [0 maxIntensity];
colormap(ax1, 'gray')

% Plot PCA image over bottom axis
ax2 = axes();
ax2.Position = ax1.Position;
plotData = combData;
plotData(end, 1:7) = 1:7;
imagesc(plotData, 'alphadata', ~isnan(combData));
axis image
%colormap(gca, [rgb('magenta'); rgb('lime'); rgb('white')])
colormap(gca, [rgb(colors{1}); rgb(colors{2}); rgb('white'); rgb(colors{3}); rgb('white'); rgb('white'); ...
        rgb('white')])
axis off



