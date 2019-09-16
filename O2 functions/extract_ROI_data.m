function extract_ROI_data(parentDir, sessionDataFile, ROIfile)
try
addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

% Load ROI info
load(fullfile(parentDir, ROIfile)); % Contains variable 'ROImetadata'
nROIs = numel(ROImetadata);
disp('ROIs loaded')
write_to_log('ROIs loaded', mfilename);

% Load imaging data
if exist(fullfile(parentDir, 'analysisMetadata.mat'), 'file')
    load(fullfile(parentDir, sessionDataFile)); % contains variable 'wholeSession'
    load(fullfile(parentDir, 'analysisMetadata.mat'));
else
    [analysisMetadata, wholeSession] = load_imaging_data(parentDir, sessionDataFile, ...
        'LoadSessionData', 1);
end

nVolumes = analysisMetadata.nVolumes;
nTrials = analysisMetadata.nTrials;
nPlanes = analysisMetadata.nPlanes;

% Extract and save xy-averaged data from each ROI
ROIDataAvg = [];
ROIDataSum = zeros(nROIs, nVolumes * nTrials);
ROIDataCount = zeros(nROIs, nVolumes * nTrials);
for iPlane = 1:nPlanes

    load(fullfile(parentDir, ['rigid_reg_sid_', num2str(analysisMetadata.sid), '_chan_1_plane_', num2str(iPlane), ...
        '_sessionFile.mat'])); % --> 'sessionFile', [y, x, volume]
    sz = size(sessionData);
    nPixels = sz(1) * sz(2);
    sessionData = reshape(sessionData, [nPixels, (nVolumes * nTrials)]); % --> [Pixel, volume]
    
    for iParent = 1:nROIs
        disp(['Extracting data for ROI #', num2str(iParent), ' of ', num2str(nROIs), '...'])
%         allROIData{iParent} = [];
        for iROI = 1:length(ROImetadata{iParent})
            currMask = ROImetadata{iParent}(iROI).mask(:); % --> [pixel]
            currPlane = ROImetadata{iParent}(iROI).plane;
            disp(currPlane);
            disp(ROImetadata);
            disp(ROImetadata{iParent})
            disp(ROImetadata{iParent}(iROI));
            if currPlane == iPlane
                ROIDataSum(iParent, :) = ROIDataSum(iParent, :) + sum(sessionData(currMask, :), 1); % --> [parentROI, volume]
                ROIDataCount(iParent) = ROIDataCount(iParent) + size(sessionData(currMask, :), 1);  % --> [parentROI]
            end
        end
    end%iParent
end%iPlane
for iParent = 1:nROIs
   currDataAvg = ROIDataSum(iParent, :) ./ ROIDataCount(iParent);       % --> [volume]
end

% CALCULATE MEAN dF/F WITHIN ROIs THROUGHOUT ENTIRE EXPERIMENT

% Using bottom 5% of entire ROI's mean value throughout each trial as baseline
ROIDataAvgSorted = sort(ROIDataAvg, 1);                                             % --> [volume, trial, ROI]
baselineMean = mean(ROIDataAvgSorted(1:round(nVolumes * 0.05), :, :), 1);           % --> [trial, ROI]
baselineMeanRep = baselineMean(ones(1, nVolumes), :, :);                            % --> [volume, trial, ROI]
ROIDffAvg = (ROIDataAvg - baselineMeanRep) ./ baselineMeanRep;                      % --> [volume, trial, ROI]

% Calculate raw fluorescence with a volume-averaged basline ROI subtracted
if nROIs > 1
    baseROIAvg = squeeze(mean(squeeze(ROIDataAvg(:,:, end)), 1));                           % --> [trial]
    ROIDataBaseSub = ROIDataAvg(:,:,end - 1) - repmat(baseROIAvg, nVolumes, 1, nROIs - 1);  % --> [volume, trial, ROI]
else
    ROIDataBaseSub = [];
end

save(fullfile(parentDir, 'ROI_Data_Avg.mat'), 'ROIDataAvg', 'ROIDffAvg', 'ROIDataBaseSub', '-v7.3') % --> [volume, trial, ROI]
catch ME
    write_to_log(getReport(ME), mfilename);
end%try
end