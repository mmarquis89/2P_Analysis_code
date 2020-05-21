saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
expList = load_expList('groupname', 'singleTrialAcq_preROIs');

expNum = 9;

% ----- AUTOMATIC PROCESSING AND CONVERSION -----
try
    
% Convert expID to old-format expDir name
currExpID = expList.expID{expNum};
disp(currExpID);
expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
expParentDir = find_parent_dir(currExpID);

% % Find experiment directory
% if isempty(dir(fullfile(parentDir, expDirName)))
%     oldExpParentDirs = dir(fullfile(parentDir, '2018 *'));
%     for iDir = 1:numel(oldExpParentDirs)
%         expParentDir = fullfile(oldExpParentDirs(iDir).folder, oldExpParentDirs(iDir).name);
%         if ~isempty(dir(fullfile(expParentDir, expDirName)))
%             break
%         end
%     end
% else
%     expParentDir = fullfile(parentDir);
% end
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
load(fullfile(expDir, 'Annotations.mat'), 'goodTrials'); 
expMd.nTrials = numel(goodTrials) + 3;

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
if size(ftDataRaw, 2) == 25
   ftDataRaw = ftDataRaw(:, 1:24);
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
for iTrial = 4:expMd.nTrials
    
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
for iTrial = 4:expMd.nTrials
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


%% ----- Save data in group analysis directory -----

try
writetable(expMd, fullfile(saveDir, [expID, '_expMetadata.csv']));
writetable(trialMd, fullfile(saveDir, [expID, '_trialMetadata.csv']));
writetable(daqOutputData, fullfile(saveDir, [expID, '_daqOutputData.csv']));
save(fullfile(saveDir, [expID, '_ficTracData.mat']), 'ftData'); 
% save(fullfile(saveDir, [expID, '_refImages.mat']), 'refImages');
save(fullfile(saveDir, [expID, '_fullExpRefImages.mat']), 'fullExpRefImages');
save(fullfile(saveDir, [expID, '_roiData.mat']), 'roiData');

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

odorANames = {'Acetoin'};
odorAConcentrations = {'e-1'};
odorBNames = {'ACV'};
odorBConcentrations = {'e-1'};
flowRates = {30};
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

disp('Odor event data saved');

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










