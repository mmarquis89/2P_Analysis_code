saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
% parentDir = 'F:\ImagingData'; 
% expList = readtable(fullfile(saveDir, 'oldExpList_singleTrialAcq_pre_ROIs.csv'), 'delimiter', ',');

expList = table({'20180623-3'}, 'variablenames', {'expID'}); 
expList.expName{1} = 'D-ANT_6s';

expNum = 1;

% ----- AUTOMATIC PROCESSING AND CONVERSION -----
try
    
% Convert expID to old-format expDir name
currExpID = expList.expID{expNum};
disp(currExpID);
expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];

% Find experiment directory
expParentDir = find_parent_dir(currExpID);
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
% expMd.nTrials = aD.nTrials;
load(fullfile(expDir, 'Annotations.mat'), 'goodTrials'); 
expMd.nTrials = numel(goodTrials);

% ----- CREATE TRIAL METADATA TABLE -----
trialMd = [];
for iTrial = 1:expMd.nTrials
    newRow = table({expID}, 'VariableNames', {'expID'});
    newRow.trialNum = iTrial;
    newRow.trialDuration = aD.trialDuration;
    newRow.nVolumes = aD.nVolumes;
    newRow.nDaqSamples = nan;
    newRow.nPanelsFrames = expMd.panelsDisplayRate * newRow.trialDuration;
    newRow.usingOptoStim = 0;
    newRow.usingPanels = 0;
    newRow.using2P = 1;
    newRow.originalTrialCount = 1;
    newRow.pmtShutoffVols = {[]};
    trialMd = [trialMd; newRow];
end

% ----- DAQ OUTPUT DATA (DOWNSAMPLED FURTHER TO 100 Hz) -----
mdFiles = dir(fullfile(expDir, 'metadata*.mat'));
if numel(mdFiles) ~= expMd.nTrials
   error('metadata nTrials mismatch'); 
end
trialNum = [];
optoStimCommand = [];
speakerCommand = [];
odorACommand = [];
odorBCommand = [];
for iFile = 1:numel(mdFiles)
   currFile = fullfile(expDir, mdFiles(iFile).name);
   disp(mdFiles(iFile).name);
   load(currFile, 'metaData');
   currOutputData = metaData.outputData;
   if iFile == 1
       expMd.daqSampRate = size(currOutputData, 1) / trialMd.trialDuration(1);
       trialNum = zeros(expMd.daqSampRate * sum(trialMd.trialDuration), 1);
       optoStimCommand = zeros(size(trialNum));
       speakerCommand = zeros(size(trialNum));
       odorACommand = zeros(size(trialNum));
       odorBCommand = zeros(size(trialNum));
       endSample = 0;
   end
   trialMd.nDaqSamples(iFile) = size(currOutputData, 1);
   
   % Column 2 = speaker, Column 3 = odorA, column 4 = odorB
   startSample = endSample + 1;
   endSample = endSample + size(currOutputData, 1);
   trialNum(startSample:endSample) = ones(size(currOutputData, 1), 1) * trialMd.trialNum(iFile);
   optoStimCommand(startSample:endSample) = zeros(size(currOutputData, 1), 1);
   speakerCommand(startSample:endSample) = currOutputData(:, 2);
   odorACommand(startSample:endSample) = currOutputData(:, 3);
   odorBCommand(startSample:endSample) = currOutputData(:, 4);
end
% Downsample
dsFactor = expMd.daqSampRate / 100;
dsInds = 1:dsFactor:size(trialNum, 1);
daqOutputData = table(trialNum(dsInds), optoStimCommand(dsInds), speakerCommand(dsInds), ...
        odorACommand(dsInds), odorBCommand(dsInds), optoStimCommand(dsInds), ...
        'VariableNames', {'trialNum', 'optoStim', 'speaker', 'odorA', 'odorB', 'optoStimSmoothed'});
expMd.daqSampRate = 100;
trialMd.nDaqSamples = trialMd.nDaqSamples / 100;
disp(currExpID);

% ----- CREATE FICTRAC DATA STRUCTURE -----

% Load .csv file with FicTrac data for all goodTrials
opts = detectImportOptions(fullfile(expDir, 'allTrials.csv'));
opts.VariableNames = {'trialNum', 'frameCounter', ...
        'dRotCamX', 'dRotCamY', 'dRotCamZ', 'dRotError', 'dRotLabX', 'dRotLabY', 'dRotLabZ', ...
        'absOrientCamX', 'absOrientCamY', 'absOrientCamZ', 'absOrientLabX', 'absOrientLabY', ...
        'absOrientCamZ', 'intX', 'intY', 'intHD', 'moveDirLab', 'moveSpeed', 'intFwMove', ...
        'intSideMove', 'timestamp', 'seqNum', 'mysteryVar1'};
ftDataRaw = readtable(fullfile(expDir, 'allTrials.csv'), opts); 
if size(ftDataRaw, 2) == 26
   ftDataRaw = ftDataRaw(:, 1:25);
end

% Check to make sure columns are in the expected positions using field with a fixed range
if abs(max(ftDataRaw.intHD) - (2*pi)) > 0.5
    error('column label mismatch')
end

% Load normalized optic flow data
flowFile = fullfile(expDir, [expDir(end-4:end), '_flow_data_norm.mat']);
if exist(flowFile, 'file')
    load(flowFile, 'flyFlowNorm');
end

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
    ftData(iTrial).yawSpeed = [0*currTrialData(1); diff(smoothdata(unwrap(ftData(iTrial).intHD), 1, ...
            'gaussian', 5, 'includenan'))];
    ftData(iTrial).fwSpeed = [0*currTrialData(1); diff(smoothdata(ftData(iTrial).intFwMove, 1, ...
            'gaussian', 5, 'includenan'))];
    ftData(iTrial).sideSpeed = [0*currTrialData(1); diff(smoothdata(ftData(iTrial).intSideMove, ...
            1, 'gaussian', 5, 'includenan'))];
    
    % Other fields
    ftData(iTrial).frameTimes = (1:aD.nFrames) / FRAME_RATE;
    ftData(iTrial).badVidFrames = isnan(ftData(iTrial).intX);
    
    % Add optic flow data
    if goodTrial && exist('flyFlowNorm', 'var')
        ftData(iTrial).meanFlow = flyFlowNorm{iTrial};
    else
        ftData(iTrial).meanFlow = nan(aD.nFrames, 1);
    end
    
end

% ----- BEHAVIOR ANNOTATION DATA -----

% Load annotation data file 
load(fullfile(expDir, 'Annotations.mat'), 'trialAnnotations'); 

quiescenceEvents = behaviorEvent('Quiescence');
isoMoveEvents = behaviorEvent('IsolatedMovement');
groomEvents = behaviorEvent('Grooming');
locEvents = behaviorEvent('Locomotion');
for iTrial = 1:expMd.nTrials
    if ~isnan(ftData(iTrial).intX(1))
        currAnnotData = trialAnnotations{iTrial};
        frameTimes = ftData(iTrial).frameTimes;
        
        % Separate according to behavior type
        qFrames = currAnnotData.actionNums == 0;
        iFrames = currAnnotData.actionNums == 4;
        gFrames = currAnnotData.actionNums == 3;
        lFrames = currAnnotData.actionNums == 2;
        
        % Append data for each trial
        if sum(qFrames) > 0
            quiescenceEvents = quiescenceEvents.append_annotation_data(expID, iTrial, qFrames, ...
                    frameTimes);
        end
        if sum(iFrames) > 0
            isoMoveEvents = isoMoveEvents.append_annotation_data(expID, iTrial, iFrames, frameTimes);
        end
        if sum(gFrames) > 0
            groomEvents = groomEvents.append_annotation_data(expID, iTrial, gFrames, frameTimes);
        end
        if sum(lFrames) > 0
            locEvents = locEvents.append_annotation_data(expID, iTrial, lFrames, frameTimes);
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
% % Pre 6/1/18
% roiData = {};
% load(fullfile(expDir, 'ROI_metadata.mat'), 'ROIdata');
% ROImetadata = ROIdata;
% for iROI = 1:numel(ROImetadata)
%     figure(1);clf;
%     nPlots = numel(ROImetadata);
%     subplotDims = numSubplots(nPlots);
%     for iPlot = 1:nPlots
%        subaxis(subplotDims(1), subplotDims(2), iPlot);
%        imshow(ROImetadata(iPlot).refImg, [0 MAX_INTENS]);
%        hold on
%        plot(ROImetadata(iPlot).xi, ROImetadata(iPlot).yi, 'linewidth', 2, 'color', 'r')
%        title(num2str(ROImetadata(iPlot).plane))
%     end
% end
% disp(['nROIs = ', num2str(numel(ROImetadata))])

% Post 6/1/18
roiData = {};
load(fullfile(expDir, 'ROI_metadata.mat'), 'ROImetadata');
for iROI = 1:numel(ROImetadata)
    figure(iROI);clf;
    currROI = ROImetadata{iROI};
    nPlots = numel(currROI);
    subplotDims = numSubplots(nPlots);
    for iPlot = 1:nPlots
       subaxis(subplotDims(1), subplotDims(2), iPlot);
       imshow(currROI(iPlot).refImg, [0 MAX_INTENS]);
       hold on
       plot(currROI(iPlot).xi, currROI(iPlot).yi, 'linewidth', 2, 'color', 'r')
    end
end
disp(['nROIs = ', num2str(numel(ROImetadata))])

% Activate figures in reverse order so first one is on top
for iROI = numel(ROImetadata):-1:1
   figure(iROI); 
end

catch ME; rethrow(ME); end

%% Reformat ROI data to match newer experiments

roiNames = {'TypeD-R', 'ANT-R'};

try
    
roiData = {};
load(fullfile(expDir, 'ROI_metadata.mat'), 'ROImetadata');
% ROImetadata = ROIdata;

load(fullfile(expDir, 'ROI_data_avg.mat'), 'ROIDataAvg'); % --> [volume, trial, ROI]
nROIs = size(ROIDataAvg, 3);

% Get raw average fluorescence data for each ROI and calculate trial-based dF/F
for iTrial = 1:expMd.nTrials
    currROIData = permute(squeeze(ROIDataAvg(:, iTrial, :)), [2 1]); % --> [ROI, volume]    
    roiData{iTrial} = struct();
    for iROI = 1:nROIs
%         currROI = ROImetadata(iROI); % pre 6/1/18
        currROI = ROImetadata{iROI}; % post 6/1/18
        
%         % Copy metadata % pre 6/1/18
%         roiData{iTrial}(iROI).name = roiNames{iROI};
%         roiData{iTrial}(iROI).subROIs = struct();
%         roiData{iTrial}(iROI).subROIs.plane = currROI.plane;
%         roiData{iTrial}(iROI).subROIs.position = [currROI.xi, currROI.yi];
% % %         
        % Copy metadata % post 6/1/18
        roiData{iTrial}(iROI).name = roiNames{iROI};
        roiData{iTrial}(iROI).subROIs = struct();
        for iSubROI = 1:numel(currROI)
            roiData{iTrial}(iROI).subROIs(iSubROI).plane = currROI(iSubROI).plane;
            roiData{iTrial}(iROI).subROIs(iSubROI).position = [currROI(iSubROI).xi, currROI(iSubROI).yi];
        end
        
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
save(fullfile(saveDir, [expID, '_trialMetadata.mat']), 'trialMd');
% writetable(trialMd, fullfile(saveDir, [expID, '_trialMetadata.csv']));
writetable(daqOutputData, fullfile(saveDir, [expID, '_daqOutputData.csv']));
save(fullfile(saveDir, [expID, '_ficTracData.mat']), 'ftData'); 
% save(fullfile(saveDir, [expID, '_refImages.mat']), 'refImages');
save(fullfile(saveDir, [expID, '_fullExpRefImages.mat']), 'fullExpRefImages');
% save(fullfile(saveDir, [expID, '_roiData.mat']), 'roiData');

if ~isempty(quiescenceEvents.eventData)
    quiescenceEvents.export_csv(saveDir, 'fileNamePrefix', expID);
end
if ~isempty(isoMoveEvents.eventData)
    isoMoveEvents.export_csv(saveDir, 'fileNamePrefix', expID);
end
if ~isempty(groomEvents.eventData)
    groomEvents.export_csv(saveDir, 'fileNamePrefix', expID);
end
if ~isempty(locEvents.eventData)
    locEvents.export_csv(saveDir, 'fileNamePrefix', expID);
end

if exist(fullfile(expDir, 'sidTrialCounts.csv'), 'file')
    copyfile(fullfile(expDir, 'sidTrialCounts.csv'), ...
            fullfile(saveDir, [expID, '_sidTrialCounts.csv'])); 
end
catch ME; rethrow(ME); end



%% ----- ODOR EVENT DATA -----

odorANames = {'EtOH'};
odorAConcentrations = {'neat'};
odorBNames = {'CO2'};
odorBConcentrations = {'4%'};
flowRates = {20};
trialNums = {[]};

try 
    
if any(daqOutputData.odorA + daqOutputData.odorB)
    odorEvents = odorEvent();
    
    for iCond = 1:numel(trialNums)
        
        % Get info for current set of conditions
        currNameA = odorANames{iCond};
        currNameB = odorBNames{iCond};
        concA = odorAConcentrations{iCond};
        concB = odorBConcentrations{iCond};
        if numel(flowRates) == 1
            currFlowRate = flowRates{1};
        else
            currFlowRate = flowRates{iCond};
        end
        if isempty(trialNums{iCond})
            currTrialNums = 1:expMd.nTrials;
        else
            currTrialNums = trialNums{iCond};
        end
        mdFieldNames = {'odorName', 'concentration', 'flowRate'};
        
        for iTrial = 1:numel(currTrialNums)
            
            currTrialData = daqOutputData(daqOutputData.trialNum == currTrialNums(iTrial), :);
            daqSampleTimes = (1:size(currTrialData, 1)) / expMd.daqSampRate;
            
            % Odor A 
            if any(currTrialData.odorA)
                
                mdFieldVals = {currNameA, concA, currFlowRate}; 
                
                % Identify onset and offset samples
                onsetSamples = find(diff(currTrialData.odorA) == 1) + 1;
                if currTrialData.odorA(1) == 1
                   onsetSamples = [1; onsetSamples]; 
                end
                offsetSamples = find(diff(currTrialData.odorA) == -1) + 1;
                if currTrialData.odorA(end) == 1
                   offsetSamples = [offsetSamples; numel(currTrialData.odorA)]; 
                end
                
                % Convert to times
                onsetTimes = daqSampleTimes(onsetSamples);
                offsetTimes = daqSampleTimes(offsetSamples);
                eventTimes = {[onsetTimes', offsetTimes']};              
                
                odorEvents = odorEvents.append_data(expMd.expID{:}, currTrialNums(iTrial), ...
                        eventTimes, mdFieldNames, mdFieldVals);
            end%if
            
            % Odor B 
            if any(currTrialData.odorB)
                
                mdFieldVals = {currNameB, concB, currFlowRate}; 
                
                % Identify onset and offset samples
                onsetSamples = find(diff(currTrialData.odorB) == 1) + 1;
                if currTrialData.odorB(1) == 1
                   onsetSamples = [1, onsetSamples]; 
                end
                offsetSamples = find(diff(currTrialData.odorB) == -1) + 1;
                if currTrialData.odorB(end) == 1
                   offsetSamples = [offsetSamples, numel(currTrialData.odorB)]; 
                end
                
                % Convert to times
                onsetTimes = daqSampleTimes(onsetSamples);
                offsetTimes = daqSampleTimes(offsetSamples);
                eventTimes = {[onsetTimes', offsetTimes']};              
                
                odorEvents = odorEvents.append_data(expMd.expID{:}, currTrialNums(iTrial), ...
                        eventTimes, mdFieldNames, mdFieldVals);
            end%if

        end%iTrial
    end%iCond
    
    % Export to .csv file
    odorEvents.export_csv(saveDir, 'fileNamePrefix', expMd.expID{:});
    
end%if

catch ME; rethrow(ME); end



%% ----- SPEAKER EVENT DATA -----

toneFreqs = {200};
ampSettings = {-30};
trialNums = {[]};

try 
    
if any(daqOutputData.speaker)
    soundEvents = soundStimEvent();
    
    for iCond = 1:numel(trialNums)
        
        % Get info for current set of conditions
        if numel(toneFreqs) == 1
            currFreq = toneFreqs{1};
        else
            currFreq = toneFreqs{iCond};
        end
        currAmpSetting = ampSettings{iCond};
        if isempty(trialNums{iCond})
            currTrialNums = 1:expMd.nTrials;
        else
            currTrialNums = trialNums{iCond};
        end
        mdFieldNames = {'toneFreq', 'ampSetting'};
        
        for iTrial = 1:numel(currTrialNums)
            
            currTrialData = daqOutputData(daqOutputData.trialNum == currTrialNums(iTrial), :);
            daqSampleTimes = (1:size(currTrialData, 1)) / expMd.daqSampRate;
            
            if any(currTrialData.speaker)
                
                mdFieldVals = {currFreq, currAmpSetting}; 
                
                % Identify onset and offset samples
                onsetSamples = find(diff(logical(currTrialData.speaker)) == 1) + 1;
                if logical(currTrialData.speaker(1)) == 1
                   onsetSamples = [1; onsetSamples]; 
                end
                offsetSamples = find(diff(logical(currTrialData.speaker)) == -1) + 1;
                if logical(currTrialData.speaker(end)) == 1
                   offsetSamples = [offsetSamples; numel(currTrialData.speaker)]; 
                end
                
                % Convert to times
                onsetTimes = daqSampleTimes(onsetSamples);
                offsetTimes = daqSampleTimes(offsetSamples);
                eventTimes = {[onsetTimes', offsetTimes']};              
                
                soundEvents = soundEvents.append_data(expMd.expID{:}, currTrialNums(iTrial), ...
                        eventTimes, mdFieldNames, mdFieldVals);
            end%if
        end%iTrial
    end%iCond
    
    % Export to .csv file
    soundEvents.export_csv(saveDir, 'fileNamePrefix', expMd.expID{:});
    
end%if

catch ME; rethrow(ME); end










