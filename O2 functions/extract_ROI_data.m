function extract_ROI_data(parentDir, sessionDataFile, ROIfile)
try
addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

% Load ROI info
load(fullfile(parentDir, ROIfile)); % Contains variable 'ROImetadata'
nROIs = numel(ROImetadata);
disp('ROIs loaded')

% Load imaging data
[analysisMetadata, wholeSession] = load_imaging_data(parentDir, sessionDataFile, 'LoadSessionData', 1);

% Extract and save xy-averaged data from each ROI
ROIDataAvg = [];
if iscell(ROImetadata)
    for iParent = 1:nROIs
        disp(['Extracting data for ROI #', num2str(iParent), ' of ', num2str(nROIs), '...'])
        currParentROIData = [];
        for iROI = 1:length(ROImetadata{iParent})
            currMask = ROImetadata{iParent}(iROI).mask;
            currPlane = ROImetadata{iParent}(iROI).plane;
            disp(size(wholeSession))
            disp(currPlane);
            disp(ROImetadata);
            disp(ROImetadata{iParent})
            disp(ROImetadata{iParent}(iROI));
            currPlaneData = squeeze(wholeSession(:,:,currPlane,:,:)); % --> [y, x, volume, trial]
            currPlaneData(~currMask(:,:,ones(1, analysisMetadata.nVolumes), ones(1, analysisMetadata.nTrials))) = nan;                    % --> [y, x, volume, trial]
            currDataLin = reshape(currPlaneData, size(currPlaneData, 1)*size(currPlaneData, 2), ...
                analysisMetadata.nVolumes, analysisMetadata.nTrials);                                   % --> [pixel, volume, trial, ROI]
            currParentROIData = cat(1, currParentROIData, currDataLin);                                 % --> [pixel, volume, trial, parentROI]
        end
        ROIDataAvg(:,:,iParent) = squeeze(mean(currParentROIData, 1, 'omitnan'));                       % --> [volume, trial, parentROI]
    end
else
    % This section for backwards compatibility
    for iROI = 1:nROIs
        
        disp(['Extracting data for ROI #', num2str(iROI), ' of ', num2str(nROIs), '...'])
        currMask = ROImetadata(iROI).mask;
        currPlane = ROImetadata(iROI).plane;
        currPlaneData = squeeze(wholeSession(:,:,currPlane,:,:));                                                   % --> [y, x, volume, trial]
        currPlaneData(~currMask(:,:,ones(1, analysisMetadata.nVolumes), ones(1, analysisMetadata.nTrials))) = nan;  % --> [y, x, volume, trial]
        currDataLin = reshape(currPlaneData, size(currPlaneData, 1)*size(currPlaneData, 2), ...
            analysisMetadata.nVolumes, analysisMetadata.nTrials);                                                   % --> [pixel, volume, trial, ROI]
        ROIDataAvg(:,:,iROI) = squeeze(mean(currDataLin, 1, 'omitnan'));                                            % --> [volume, trial, ROI]
    end
end

% CALCULATE MEAN dF/F WITHIN ROIs THROUGHOUT ENTIRE EXPERIMENT

% Using bottom 5% of entire ROI's mean value throughout each trial as baseline
ROIDataAvgSorted = sort(ROIDataAvg, 1);                                                      % --> [volume, trial, ROI] 
baselineMean = mean(ROIDataAvgSorted(1:round(analysisMetadata.nVolumes * 0.05), :, :), 1);   % --> [trial, ROI] 
baselineMeanRep = baselineMean(ones(1, analysisMetadata.nVolumes), :, :);                    % --> [volume, trial, ROI] 
ROIDffAvg = (ROIDataAvg - baselineMeanRep) ./ baselineMeanRep;                               % --> [volume, trial, ROI] 

save(fullfile(parentDir, 'ROI_Data_Avg.mat'), 'ROIDataAvg', 'ROIDffAvg', '-v7.3')            % --> [volume, trial, ROI]
catch ME
    write_to_log(getReport(ME), mfilename);
end%try
end