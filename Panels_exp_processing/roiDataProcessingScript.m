
parentDir = 'D:\Dropbox (HMS)\2P Data\20191107-1_38A11-Chrimson_60D05-7f\ProcessedData';
roiDefFileStr = 'roiDefs_*trial*.mat';
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
    
    disp(allROIData(iFile).trialNum);
    
    % Load current trial's imaging data
    currImgDataFile = imgDataFiles(imgDataTrialNums == allROIData(iFile).trialNum);
    load(fullfile(parentDir, currImgDataFile.name)); % --> 'imageData' [y, x, plane, volume]
    allROIData(iFile).imgDataFile = currImgDataFile.name;
    
% imgData = imgDataReg; clear imgDataReg;
    
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
        allROIData(iFile).roiDefs(iROI).data = mean(currRoiData, 1); % --> [volume]
        
    end%iROI
    
    clear imgData;
    
end%iFile


% Save file with all ROI data
save(fullfile(parentDir, [saveFileStr, '.mat']), 'allROIData', '-v7.3');




