function analysis_processing(expDate, sid)


addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

% Define directory paths
imgDataDir= ['/n/scratch2/mjm60/', expDate, '/ImagingData'];
vidDataDir= ['/n/scratch2/mjm60/', expDate, '/BehaviorVideo'];
imgSaveDir= ['/n/scratch2/mjm60/', expDate, '/sid_', num2str(sid), '/ImagingData'];
vidSaveDir= ['/n/scratch2/mjm60/', expDate, '/sid_', num2str(sid), '/BehaviorVideo'];

sessionDataFile = ['rigid_sid_', num2str(sid), '_sessionFile.mat'];

% Initialize cluster communication
c = parcluster; 
% 
% Save structure of analysis metadata + hardcoded parameters ('analysisMetadata.mat')
[analysisMetadata, ~] = load_imaging_data(imgSaveDir, sessionDataFile);
disp('Output:')
disp(analysisMetadata)
analysisMetadata.ROIdata = [];
analysisMetadata.MAX_INTENSITY = 1000; % To control brightness of ref image plots

% Load real-time FicTrac data, downsample and add to analysis metadata
ftDataFiles = dir(fullfile(vidDataDir, ['fictrac*sid_*', num2str(sid), '_bid*.mat']));
fileNames = sort({ftDataFiles.name})';
analysisMetadata.daqFtData = [];
for iFile = 1:numel(fileNames)
    load(fullfile(vidDataDir, fileNames{iFile})); % 'blockData' --> each col is an analog input channel
    currFt = blockData;
    analysisMetadata.blockData(iFile).inputData = currFt;
    sampPerTrial = size(currFt, 1) / analysisMetadata.blockData(iFile).nTrials;
    rsFtData = reshape(currFt', size(currFt, 2), sampPerTrial, analysisMetadata.blockData(iFile).nTrials);
    disp(size(rsFtData));
    analysisMetadata.daqFtData = cat(3, analysisMetadata.daqFtData, rsFtData(:,1:100:end,:));
end

% Process and save annotation types ('annotationTypes.mat')
save(fullfile(imgSaveDir, 'analysisMetadata.mat'), 'analysisMetadata', '-v7.3')
[annotationTypes, annotationTypeSummary] = process_annotation_types(analysisMetadata);
save(fullfile(imgSaveDir, 'annotationTypes.mat'), 'annotationTypes', 'annotationTypeSummary', '-v7.3');
write_to_log('Annotation types saved', mfilename);

m = matfile(fullfile(imgSaveDir, sessionDataFile));
sz = size(m, 'wholeSession');
clear m
nPixels = prod(sz(1:2));

nTrials = analysisMetadata.nTrials; 
nVolumes = analysisMetadata.nVolumes; 
nPlanes = analysisMetadata.nPlanes;

% Calculate overall behavior state dF/F
memGB = ceil(1.5e-08 * nPixels * nTrials * nVolumes * nPlanes);
if memGB > 249
    memGB = 249;
end
timeLimitMin = ceil(1e-08 * nPixels * nTrials * nVolumes * nPlanes);
if timeLimitMin > 719
    timeLimitMin = 719;
end
queueName = 'short';
jobName = [expDate, '_sid_', num2str(sid), '_behavior_state_dFF_calc'];
c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
inputArgs = {imgSaveDir, sessionDataFile};
behavStateDffCalcJob = c.batch(@behavioral_state_dff_calc, 0, inputArgs);

% Extract ROI data
ROIfile = 'ROI_metadata.mat';
sessionDataFile = ['rigid_sid_', num2str(sid), '_sessionFile.mat'];
memGB = ceil(2.5e-08 * nPixels * nTrials * nVolumes * nPlanes);
if memGB > 249
    memGB = 249;
end
timeLimitMin = ceil(2e-08 * nPixels * nTrials * nVolumes * nPlanes);
if timeLimitMin > 719
    timeLimitMin = 719;
end
queueName = 'short';
jobName = [expDate, '_sid_', num2str(sid), '_extract_ROI_data'];
c = set_job_params(c, queueName, timeLimitMin, memGB, jobName); 
inputArgs = {imgSaveDir, sessionDataFile, ROIfile};
c.batch(@extract_ROI_data, 0, inputArgs);


end%function