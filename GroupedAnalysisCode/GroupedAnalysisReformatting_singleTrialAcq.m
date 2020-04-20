saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2018'; 
expList = readtable(fullfile(saveDir, 'oldExpList_singleTrialAcq.csv'), 'delimiter', ',');

expNum = 20;

% ----- AUTOMATIC PROCESSING AND CONVERSION -----
try
    
% Convert expID to old-format expDir name
currExpID = expList.expID{expNum};
disp(currExpID);
expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];

% Find experiment directory
if isempty(dir(fullfile(parentDir, expDirName)))
    oldExpParentDirs = dir(fullfile(parentDir, '2018 *'));
    for iDir = 1:numel(oldExpParentDirs)
        expParentDir = fullfile(oldExpParentDirs(iDir).folder, oldExpParentDirs(iDir).name);
        if ~isempty(dir(fullfile(expParentDir, expDirName)))
            break
        end
    end
else
    expParentDir = fullfile(parentDir);
end
sidDir = dir(fullfile(expParentDir, expDirName, 'sid_*'));
sidDir = sidDir([sidDir.isdir]);
if numel(sidDir) == 1
    expDir = fullfile(expParentDir, expDirName, sidDir.name); % If there's only one sid directory
elseif numel(sidDir) > 1
    expDir = fullfile(expParentDir, expDirName, 'sid_master'); % If I've already manually combined sids
else
    error('Could not find experiment directory');
end
expName = expList.expName{expNum};

% Load analysis metadata file
load(fullfile(expDir, 'analysisMetadata.mat'), 'analysisMetadata');
aD = analysisMetadata;

FRAME_RATE = 25;% Behavior video (and therefore FicTrac) frame rate was a constant;

% ----- CREATE EXP METADATA TABLE -----
expID = regexprep(aD.expDate, {'_', 'exp'}, {'', '-'});
if ~strcmp(expID(end - 1), '-')
   expID = [expID, '-1']; % Add exp num of 1 if there wasn't one in the expDate
end
expMd = table({expID}, 'VariableNames', {'expID'});
expMd.expName = {expName};
expMd.daqSampRate = nan;
expMd.panelsDisplayRate = 50;
expMd.volumeRate = aD.volumeRate;
expMd.nPlanes = aD.nPlanes;
expMd.nTrials = aD.nTrials;


% ----- CREATE TRIAL METADATA TABLE -----
trialMd = [];
for iTrial = 1:expMd.nTrials
    newRow = table(iTrial, 'VariableNames', {'trialNum'});
    newRow.trialDuration = aD.trialDuration;
    newRow.nVolumes = aD.nVolumes;
    newRow.nDaqSamples = nan;
    newRow.nPanelsFrames = expMd.panelsDisplayRate * newRow.trialDuration;
    newRow.usingOptoStim = 0;
    newRow.usingPanels = 0;
    newRow.using2P = 1;
    trialMd = [trialMd; newRow];
end


% ----- CREATE FICTRAC DATA STRUCTURE -----

% Load .csv file with FicTrac data for all goodTrials
opts = detectImportOptions(fullfile(expDir, 'allTrials.csv'));
opts.VariableNames = {'trialNum', 'frameCounter', ...
        'dRotCamX', 'dRotCamY', 'dRotCamZ', 'dRotError', 'dRotLabX', 'dRotLabY', 'dRotLabZ', ...
        'absOrientCamX', 'absOrientCamY', 'absOrientCamZ', 'absOrientLabX', 'absOrientLabY', ...
        'absOrientCamZ', 'intX', 'intY', 'intHD', 'moveDirLab', 'moveSpeed', 'intFwMove', ...
        'intSideMove', 'timestamp', 'seqNum', 'mysteryVar'};
ftDataRaw = readtable(fullfile(expDir, 'allTrials.csv'), opts); 

% Check to make sure columns are in the expected positions using fields with a fixed range
if abs(max(ftDataRaw.absOrientCamZ) - pi) > 0.5 || abs(max(ftDataRaw.intHD) - (2*pi)) > 0.5
    error('column label mismatch')
end

% Load normalized optic flow data
load(fullfile(expDir, [expDir(end-4:end), '_flow_data_norm.mat']), 'flyFlowNorm');

% Build structure
ftData = struct();
for iTrial = 1:expMd.nTrials
    
    % Add fields from raw data
    currTrialData = ftDataRaw{ftDataRaw.trialNum == iTrial, :};
    if size(currTrialData, 1) ~= aD.nFrames
        currTrialData = nan(aD.nFrames, size(currTrialData, 2));
        currTrialData(:, 1) = iTrial;
        goodTrial = 0;
    else
        goodTrial = 1;
    end
    ftData(iTrial).intX = currTrialData(:, 16);
    ftData(iTrial).intY = currTrialData(:, 17);
    ftData(iTrial).intHD = currTrialData(:, 18);
    ftData(iTrial).moveSpeed = currTrialData(:, 20);
    ftData(iTrial).intFwMove = currTrialData(:, 21);
    ftData(iTrial).intSideMove = currTrialData(:, 22);
    
    % Calculated fields (zero multiplication is in case ftData is 'nan')
    ftData(iTrial).yawSpeed = [0*currTrialData(1), diff(smoothdata(unwrap(ftData(iTrial).intHD), 1, ...
            'gaussian', 5, 'includenan'))];
    ftData(iTrial).fwSpeed = [0*currTrialData(1), diff(smoothdata(ftData(iTrial).intFwSpeed), 1, ...
            'gaussian', 5, 'includenan')];
    ftData(iTrial).sideSpeed = [0*currTrialData(1), diff(smoothdata(ftData(iTrial).intSideSpeed), ...
            1, 'gaussian', 5, 'includenan')];
    
    % Other fields
    ftData(iTrial).frameTimes = (1:trialMd(iTrial).trialDuration) / FRAME_RATE;
    ftData(iTrial).badVidFrames = isnan(ftData(iTrial).intX);
    
    % Add optic flow data
    if goodTrial
        ftData(iTrial).meanFlow = flyFlowNorm{iTrial};
    else
        ftData(iTrial).meanFlow = nan(aD.nFrames, 1);
    end
    
end

% ----- BEHAVIOR ANNOTATION DATA -----

% Load annotation data file 
load(fullfile(expDir, 'Annotations.mat'), 'trialAnnotations'); 

quiescenceEvents = behaviorEvent(expID, 'Quiescence');
isoMoveEvents = behaviorEvent(expID, 'IsolatedMovement');
groomEvents = behaviorEvent(expID, 'Grooming');
locEvents = behaviorEvent(expID, 'Locomotion');
for iTrial = 1:expMd.nTrials
   if ~isnan(ftData(iTrial).intX(1))
       currAnnotData = trialAnnotations{iTrial};
       frameTimes = ftData(iTrial).frameTimes;
       
       % Separate according to behavior type
       qFrames = currAnnotData == 0;
       iFrames = currAnnotData == 4;
       gFrames = currAnnotData == 3;
       lFrames = currAnnotData == 2;
       
       % Append data for each trial
       if sum(qFrames) > 0
           quiescenceEvents = quiescenceEvents.append_annotation_data(iTrial, qFrames, frameTimes);
       end
       if sum(iFrames) > 0
           isoMoveEvents = isoMoveEvents.append_annotation_data(iTrial, iFrames, frameTimes);
       end
       if sum(gFrames) > 0
           groomEvents = groomEvents.append_annotation_data(iTrial, gFrames, frameTimes);
       end
       if sum(lFrames) > 0
           locEvents = locEvents.append_annotation_data(iTrial, lFrames, frameTimes);
       end
   end
end%iTrial


% ----- REFERENCE IMAGES -----

% Just using the full experiment reference images for all the trials
refImages = [];
for iPlane = 1:expMd.nPlanes
    fullExpRefImages(:, :, iPlane) = aD.refImg{iPlane};
end

catch ME; rethrow(ME); end

%% Generate ROI names

MAX_INTENS = 700;

try
    
roiData = {};
load(fullfile(expDir, 'ROI_metadata.mat'), 'ROImetadata');
for iROI = 1:numel(ROImetadata)
    figure(1);clf;
    nPlots = numel(ROImetadata);
    subplotDims = numSubplots(nPlots);
    for iPlot = 1:nPlots
       subaxis(subplotDims(1), subplotDims(2), iPlot);
       imshow(ROImetadata(iPlot).refImg, [0 MAX_INTENS]);
       hold on
       plot(ROImetadata(iPlot).xi, ROImetadata(iPlot).yi, 'linewidth', 2, 'color', 'r')
    end
end
disp(['nROIs = ', num2str(numel(ROImetadata))])

catch ME; rethrow(ME); end

%% Reformat ROI data to match newer experiments

roiNames = {'TypeD-L', 'ANT-L', 'LH-L'};

try
    
roiData = {};
load(fullfile(expDir, 'ROI_metadata.mat'), 'ROImetadata');

load(fullfile(expDir, 'ROI_data_avg.mat'), 'ROIDataAvg'); % --> [volume, trial, ROI]
nROIs = size(ROIDataAvg, 3);

% Get raw average fluorescence data for each ROI and calculate trial-based dF/F
for iTrial = 1:expMd.nTrials
    currROIData = permute(squeeze(ROIDataAvg(:, iTrial, :)), [2 1]); % --> [ROI, volume]    
    roiData{iTrial} = struct();
    for iROI = 1:nROIs
        currROI = ROImetadata(iROI);
        
        % Copy metadata
        roiData{iTrial}(iROI).name = roiNames{iROI};
        roiData{iTrial}(iROI).subROIs = struct();
        roiData{iTrial}(iROI).subROIs.plane = currROI.plane;
        roiData{iTrial}(iROI).subROIs.position = [currROI.xi, currROI.yi];
        
        % Copy raw avg fluorescence for current ROI
        roiData{iTrial}(iROI).rawFl = currROIData(iROI, :);
        
        % Calculate trial-based dF/F
        roiDataSorted = sort(roiData{iTrial}(iROI).rawFl);
        baselineF = median(roiDataSorted(1:round(numel(roiDataSorted) * 0.05))); % Bottom 5% as baseline
        roiData{iTrial}(iROI).dffData = (roiData{iTrial}(iROI).rawFl - baselineF) ./ baselineF;        
        
    end
end

% Get list of all unique ROI names in experiment
allROINames = {};
for iTrial = 1:numel(roiData)
    allROINames = [allROINames, {roiData{iTrial}.name}];
end
ROIList = unique(allROINames);

% Extract and concatenate all fl data throughout experiment for each ROI
rawROIData = repmat({[]}, 1, numel(ROIList)); stdDevs = [];
for iROI = 1:numel(ROIList)
    for iTrial = 1:numel(roiData)
        if sum(strcmp({roiData{iTrial}.name}, ROIList{iROI}))
            
            currROIData = roiData{iTrial}(strcmp({roiData{iTrial}.name}, ...
                    ROIList{iROI})).rawFl;
            
            FL_THRESH = 8; % Hack to prevent the inclusion of trials when the PMT shutter was closed
            stdDevs(iROI, iTrial) = std(currROIData);
            if std(currROIData) > FL_THRESH
                rawROIData{iROI} = [rawROIData{iROI}, currROIData];
            end
            
        end
    end%iTrial
end%iROI

% Calculate a whole-experiment baseline F value for each ROI
for iROI = 1:numel(ROIList)
    currFlSorted = sort(rawROIData{iROI});
    currROIBaseline = median(currFlSorted(1:round(numel(currFlSorted) * 0.05))); % This is actually just the 2.5-th percentile I see?
    
    % Calculate experiment-wide dF/F for each trial
    for iTrial = 1:numel(roiData)
         if sum(strcmp({roiData{iTrial}.name}, ROIList{iROI}))
             expDffData = (roiData{iTrial}(strcmp({roiData{iTrial}.name}, ...
                    ROIList{iROI})).rawFl - currROIBaseline)./ currROIBaseline;
             roiData{iTrial}(strcmp({roiData{iTrial}.name}, ...
                    ROIList{iROI})).expDffData = expDffData;
         end
    end    
end

catch ME; rethrow(ME); end

%% ----- Save data in group analysis directory -----
try
writetable(expMd, fullfile(saveDir, [expID, '_expMetadata.csv']));
writetable(trialMd, fullfile(saveDir, [expID, '_trialMetadata.csv']));
% writetable(daqOutputData, fullfile(saveDir, [expID, '_daqOutputData.csv']));
save(fullfile(saveDir, [expID, '_ficTracData.mat']), 'ftData'); 
% save(fullfile(saveDir, [expID, '_refImages.mat']), 'refImages');
save(fullfile(saveDir, [expID, '_fullExpRefImages.mat']), 'fullExpRefImages');
save(fullfile(saveDir, [expID, '_roiData.mat']), 'roiData');

if ~isempty(quiescenceEvents.eventData)
    quiescenceEvents.export_csv(saveDir, '');
end
if ~isempty(isoMoveEvents.eventData)
    isoMoveEvents.export_csv(saveDir, '');
end
if ~isempty(groomEvents.eventData)
    groomEvents.export_csv(saveDir, '');
end
if ~isempty(locEvents.eventData)
    locEvents.export_csv(saveDir, '');
end

if exist(fullfile(expDir, 'sidTrialCounts.csv'), 'file')
    copyfile(fullfile(expDir, 'sidTrialCounts.csv'), ...
            fullfile(saveDir, [expID, '_sidTrialCounts.csv'])); 
end
catch ME; rethrow(ME); end























