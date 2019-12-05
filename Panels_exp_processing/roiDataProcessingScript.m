
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\20191203-1_38A11_ChR_60D05_7f\ProcessedData';
roiDefFileStr = 'roiDefs_trial*.mat';
imgFileStr = 'imagingData_reg*trial*.mat';
saveFileStr = 'roiData_reg'; 

% Identify imaging data files and their trial numbers
imgDataFiles = dir(fullfile(parentDir, imgFileStr));
imgDataTrialNums = get_trialNum({imgDataFiles.name});

% Identify and load ROI def files
roiDefFiles = dir(fullfile(parentDir, roiDefFileStr));
allROIData = struct;
for iFile = 1:numel(roiDefFiles)
    load(fullfile(parentDir, roiDefFiles(iFile).name)); % 'roiDefs'
    allROIData(iFile).roiDefs = roiDefs;
    allROIData(iFile).trialNum = get_trialNum(roiDefFiles(iFile).name);
    allROIData(iFile).roiDefFile = roiDefFiles(iFile).name;
end

% Extract ROI data
for iFile = 1:numel(allROIData)
    
    disp(['Processing ROI #', num2str(allROIData(iFile).trialNum), '...']);
    
    % Load current trial's imaging data
    currImgDataFile = imgDataFiles(imgDataTrialNums == allROIData(iFile).trialNum);
    load(fullfile(parentDir, currImgDataFile.name)); % --> 'imageData' [y, x, plane, volume]
    allROIData(iFile).imgDataFile = currImgDataFile.name;
        
    % Reshape into 1D frames
    sz = size(imgData);
    nVolumes = sz(4);
    nPx = sz(1) * sz(2);
    imgData = reshape(imgData, [nPx, sz(3), sz(4)]); % --> [pixel, plane, volume]
    
    % Loop through and extract data for each subROI
    for iROI = 1:numel(allROIData(iFile).roiDefs)
        currROI = allROIData(iFile).roiDefs(iROI);
        currRoiData = [];
        for iSubROI = 1:numel(currROI.subROIs)
            currSubROI = currROI.subROIs(iSubROI);
            currImgData = squeeze(imgData( :, currSubROI.plane, :));% --> [pixel, volume]
            mask = currSubROI.mask(:);                              % --> [pixel]
            currImgData = currImgData(mask, :);                     % --> [pixel, volume]
            currRoiData = [currRoiData; currImgData];               % --> [pixel, volume]
        end
        
        % Average data across pixels and subROIs
        roiDataAvg = mean(currRoiData, 1);
        allROIData(iFile).roiDefs(iROI).rawData = roiDataAvg; % --> [volume]
        
        % Calculate and save dF/F
        roiDataSorted = sort(roiDataAvg);
        baselineF = median(roiDataSorted(1:round(numel(roiDataSorted) * 0.05))); % Bottom 5% as baseline
        allROIData(iFile).roiDefs(iROI).dffData = (roiDataAvg - baselineF) ./ baselineF;
        
        % Get z-scored data
        allROIData(iFile).roiDefs(iROI).zscoreData = zscore((roiDataAvg - baselineF) ./ baselineF);
        
    end%iROI
    
    clear imgData;
    
end%iFile


% Save file with all ROI data
disp('Saving data from all ROIs...')
save(fullfile(parentDir, [saveFileStr, '.mat']), 'allROIData', '-v7.3');
disp('Saving complete!')



