
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';

% Find all expMetadata files and concatenate into a table
expMdFiles = dir(fullfile(parentDir, '*_expMetadata.csv'));
gaplessAcqExpMd = [];
for iExp = 1:numel(expMdFiles)
    currExpMd = readtable(fullfile(parentDir, expMdFiles(iExp).name));
    gaplessAcqExpMd = [gaplessAcqExpMd; currExpMd];
end
gaplessAcqExpList = gaplessAcqExpMd.expID;

% Load and concatenate data from all experiments in parent directory
allTrialMd = [];
allExpOdorEvents = soundStimEvent();
roiDataTable = [];
ftDataTable = [];
for iExp = 1:numel(gaplessAcqExpList)
    
    currExpID = gaplessAcqExpList{iExp};
    disp(currExpID);
    
    % trialMetadata files
    currTrialMd = readtable(fullfile(parentDir, [currExpID, '_trialMetadata.csv']));
    if ~any(strcmp(fieldnames(currTrialMd), 'originalTrialCount'))
       currTrialMd.originalTrialCount = ones(size(currTrialMd, 1), 1); 
    end
    allTrialMd = [allTrialMd; currTrialMd];
    
    % Odor event data
%     odorEventFile = fullfile(parentDir, [currExpID, '_event_data_odor.csv']);
    odorEventFile = fullfile(parentDir, [currExpID, '_event_data_soundstim.csv']);
    if exist(odorEventFile, 'file')
        allExpOdorEvents = allExpOdorEvents.load_csv(parentDir, 'fileNamePrefix', currExpID);
    end
%     
end%iExp


%%
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';
currExpList = unique(allExpOdorEvents.eventData.expID);
allExpRoiData = [];
for iExp = 1:numel(currExpList)
    if exist(fullfile(parentDir, [currExpList{iExp}, '_roiData.mat']))
        load(fullfile(parentDir, [currExpList{iExp}, '_roiData.mat']), 'roiData'); 
        allExpRoiData = [allExpRoiData; roiData];
    end
end

%%
targetRoiRegex = '^ANT(-.)?';
analysisWin = [2 5];

currStimData = allExpOdorEvents.eventData;
trialList = unique(currStimData(:, {'expID', 'trialNum'}));
currRoiData = innerjoin(trialList, allExpRoiData(~cellfun(@isempty, regexp(allExpRoiData.roiName, ...
        targetRoiRegex)), :));

currMd = innerjoin(innerjoin(innerjoin(trialList, gaplessAcqExpMd(:, {'expID', 'volumeRate'})), ...
        allTrialMd(:, {'expID', 'trialNum', 'trialDuration', 'nVolumes', 'originalTrialCount'})), ...
        currRoiData);

expIDs = unique(currMd.expID);
outputTable = table(expIDs, 'VariableNames', {'expID'});
for iExp = 1:numel(expIDs)
    currExpID = expIDs{iExp, 1};
    currExpMd = currMd(strcmp(currMd.expID, currExpID), :);
    nStims = size(innerjoin(currExpMd(:, {'expID', 'trialNum'}), currStimData), 1);
    winSizeVolumes = ceil(sum(analysisWin) * currExpMd.volumeRate(1)) + 1;
    currExpDff = nan(winSizeVolumes, nStims);
    stimCount = 1;
    for iTrial = 1:size(currExpMd, 1)
        currTrialDff = currExpMd.expDffData{iTrial};
        currTrialDff = reshape(currTrialDff, currExpMd.originalTrialCount(iTrial), ...
                currExpMd.nVolumes(iTrial) / currExpMd.originalTrialCount(iTrial))';
        currTrialDff = currTrialDff(:);
        volTimes = (1:currExpMd.nVolumes(iTrial)) / currExpMd.volumeRate(1);
        currTrialOdorData = innerjoin(currExpMd(iTrial, {'expID', 'trialNum'}), currStimData);
        for iStim = 1:size(currTrialOdorData, 1)
            startTime = currTrialOdorData.onsetTime(iStim) - analysisWin(1);
            endTime = currTrialOdorData.onsetTime(iStim) + analysisWin(2);
            startFrame = argmin(abs(volTimes - startTime));
            endFrame = argmin(abs(volTimes - endTime));
%             startFrame = argmin(abs(testVolTimes - startTime));
%             endFrame = argmin(abs(testVolTimes - endTime));
            currStimDff = currTrialDff(startFrame:endFrame);
            stimStartVol = 1 + argmin(abs(volTimes - analysisWin(1)));
            currStimDur = currTrialOdorData.offsetTime(iStim) - currTrialOdorData.onsetTime(iStim);
            stimEndVol = 1 + argmin(abs(volTimes - (analysisWin(1) + currStimDur)));
            currStimDff(stimStartVol) = nan;
            currStimDff(stimEndVol) = nan;
            currExpDff(1:numel(currStimDff), stimCount) = currStimDff';
            stimCount = stimCount + 1;
        end
    end
    outputTable.expDff{iExp} = currExpDff;
end

%%

for iExp = 1:size(outputTable, 1)
   
    figure(iExp);clf;
    nanVals = isnan(outputTable.expDff{iExp});
    smData = smoothdata(outputTable.expDff{iExp}, 1, 'gaussian', 3);
    smData(nanVals) = nan;
    plot(smData);
    hold on
    plot(mean(smData, 2, 'omitnan'), 'color', 'k', 'linewidth', ....
            2);
    
    
end

%% Make list of number of experiments for each unique odor identity/concentration combo

uniqueStims = unique(allExpOdorEvents.eventData(:, 5:6));
disp(uniqueStims)

odorCounts = table(nan(size(uniqueStims, 1), 1), 'VariableNames', {'count'});
for iStim = 1:size(uniqueStims, 1)
    currStims = innerjoin(uniqueStims(iStim, :), allExpOdorEvents.eventData);        
    odorCounts{iStim, 1} = size(currStims, 1);
end

disp([uniqueStims, odorCounts])

%% Extract all ROI fluorescence data in a window around a particular stimulus

targetStim = uniqueStims(4,:);
targetRoiRegex = '^ANT(-.)?';
analysisWin = [2 2];

currStimOdorData = innerjoin(targetStim, allExpOdorEvents.eventData);
trialList = unique(currStimOdorData(:, {'expID', 'trialNum'}));
currRoiData = innerjoin(trialList, roiDataTable(~cellfun(@isempty, ...
        regexp(roiDataTable{:, {'roiName'}}, targetRoiRegex, 'match', 'once')), :));
currMd = innerjoin(innerjoin(innerjoin(trialList, gaplessAcqExpMd(:, {'expID', 'volumeRate'})), ...
        allTrialMd(:, {'expID', 'trialNum', 'trialDuration', 'nVolumes', 'originalTrialCount'})), ...
        currRoiData);

expIDs = unique(currMd.expID);
outputTable = table(expIDs, 'VariableNames', {'expID'});
for iExp = 1:numel(expIDs)
    currExpID = expIDs{iExp, 1};
    currExpMd = currMd(strcmp(currMd.expID, currExpID), :);
    nStims = size(innerjoin(currExpMd(:, {'expID', 'trialNum'}), currStimOdorData), 1);
    winSizeVolumes = ceil(sum(analysisWin) * currExpMd.volumeRate(iTrial)) + 1;
    currExpDff = nan(winSizeVolumes, nStims);
    stimCount = 1;
    for iTrial = 1%:size(currExpMd, 1)
        currTrialDff = currExpMd.expDffData{iTrial};
        currTrialDff = reshape(currTrialDff, currExpMd.originalTrialCount(iTrial), ...
                currExpMd.nVolumes(iTrial) / currExpMd.originalTrialCount(iTrial))';
        currTrialDff = currTrialDff(:);
        volTimes = (1:currExpMd.nVolumes(iTrial)) / currExpMd.volumeRate(1);
        currTrialOdorData = innerjoin(currExpMd(iTrial, {'expID', 'trialNum'}), currStimOdorData);
        for iStim = 1:size(currTrialOdorData, 1)
            startTime = currTrialOdorData.onsetTime(iStim) - analysisWin(1);
            endTime = currTrialOdorData.onsetTime(iStim) + analysisWin(2);
%             startFrame = argmin(abs(volTimes - startTime));
%             endFrame = argmin(abs(volTimes - endTime));
            startFrame = argmin(abs(testVolTimes - startTime));
            endFrame = argmin(abs(testVolTimes - endTime));
            currStimDff = currTrialDff(startFrame:endFrame);
            currExpDff(1:numel(currStimDff), stimCount) = currStimDff';
            currSti
            stimCount = stimCount + 1;
        end
    end
    outputTable.expDff{iExp} = currExpDff;
end





%% Load ROI data to compare offset across trials


% PRE-GAPLESS ACQUISITION
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\singleTrialAcq_preROIs'; 
roiFiles = dir(fullfile(parentDir, '*_roiData.mat'));

allExpROIData = [];
for iExp = 1:numel(roiFiles)
    currFileName = roiFiles(iExp).name;
    currExpID = currFileName(1:10);
    disp(currExpID)
    newRow = table({currExpID}, 'VariableNames', {'expID'});
%     % Load trial metadata
%     trialMd = readtable(fullfile(saveDir, [currExpID, '_trialMetadata.csv']), 'delimiter', ',');
    
    % Load ROI data
    load(fullfile(parentDir, currFileName), 'roiData');
    roiNames = unique(roiData.roiName);
    nROIs = numel(roiNames);
    for iROI = 1:nROIs
        newRow.roiName = roiNames(iROI);
        currROIData = roiData(strcmp(roiNames{iROI}, roiData.roiName), :);
        fullExpFlData = [];
        trialBounds = 1;
        for iTrial = 1:size(currROIData, 1)
            fullExpFlData = [fullExpFlData, currROIData.rawFl{iTrial}];
            trialBounds(end + 1) = trialBounds(end) + numel(currROIData.rawFl{iTrial});
        end
        newRow.rawFl = {fullExpFlData};
        newRow.trialBounds = {trialBounds(1:end-1)};
        allExpROIData = [allExpROIData; newRow];
    end
end


% OLD FORMAT, GAPLESS ACQUISITION
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\gaplessAcq'; 
roiFiles = dir(fullfile('D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData', '*_roiData.mat'));

allExpROIData = [];
for iExp = 1:numel(roiFiles)
    currFileName = roiFiles(iExp).name;
    currExpID = currFileName(1:10);
    disp(currExpID)
    newRow = table({currExpID}, 'VariableNames', {'expID'});
    
    % Load trial metadata
    trialMd = readtable(fullfile(parentDir, [currExpID, '_trialMetadata.csv']), 'delimiter', ',');
    
    % Load ROI data
%     load(fullfile(parentDir, currFileName), 'roiData');
    load(fullfile('D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData', currFileName), 'roiData');
    roiData = innerjoin(roiData, trialMd(:, {'expID', 'trialNum', 'originalTrialCount'}));
    roiNames = unique(roiData.roiName);
    nROIs = numel(roiNames);
    for iROI = 1:nROIs
        newRow.roiName = roiNames(iROI);
        currROIData = roiData(strcmp(roiNames{iROI}, roiData.roiName), :);
        fullExpFlData = [];
        trialBounds = 1;
        for iTrial = 1:size(currROIData, 1)
            volCount = numel(currROIData.rawFl{iTrial}) / currROIData.originalTrialCount(iTrial);
            fullExpFlData = [fullExpFlData, currROIData.rawFl{iTrial}'];
            trialBounds = [trialBounds, trialBounds(end) + ...
                    (volCount:volCount:numel(currROIData.rawFl{iTrial}))];
        end
        newRow.rawFl = {fullExpFlData};
        newRow.trialBounds = {trialBounds(1:end-1)};
        allExpROIData = [allExpROIData; newRow];
    end
end

% NEW FORMAT
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\newData'; 
roiFiles = dir(fullfile(parentDir, '*_roiData.mat'));

allExpROIData = [];
for iExp = 1:numel(roiFiles)
    currFileName = roiFiles(iExp).name;
    currExpID = currFileName(1:10);
    disp(currExpID)
    newRow = table({currExpID}, 'VariableNames', {'expID'});
    
    % Load ROI data
    load(fullfile(parentDir, currFileName), 'roiData');
    roiNames = unique(roiData.roiName);
    nROIs = numel(roiNames);
    for iROI = 1:nROIs
        newRow.roiName = roiNames(iROI);
        currROIData = roiData(strcmp(roiNames{iROI}, roiData.roiName), :);
        fullExpFlData = [];
        trialBounds = 1;
        for iTrial = 1:size(currROIData, 1)
            fullExpFlData = [fullExpFlData, currROIData.rawFl{iTrial}];
            trialBounds(end + 1) = trialBounds(end) + numel(currROIData.rawFl{iTrial});
        end
        newRow.rawFl = {fullExpFlData};
        newRow.trialBounds = {trialBounds(1:end-1)};
        allExpROIData = [allExpROIData; newRow];
    end
end

%% PLOT ROI DATA TO CHECK FOR OFFSETS AT TRIAL BOUNDS

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\testing';

currFileName = 'newFormat.mat'; %'gaplessAcq.mat'; %'preGaplessAcq.mat'; %  

smWin = 7;

load(fullfile(parentDir, currFileName), 'allExpROIData');
expIDs = unique(allExpROIData.expID);
for iExp = 1:numel(expIDs)
   currData = allExpROIData(strcmp(allExpROIData.expID, expIDs{iExp}), :); 
   
   % Create figure
   f = figure(iExp); clf; hold on
   f.Color = [1 1 1];
   f.Position = [50, 250, 1500, 700];
   
   for iROI = 1:size(currData, 1)
       flData = currData.rawFl{iROI}(:);
       flData(currData.trialBounds{1}) = nan;
       plot(smoothdata(flData, 'gaussian', smWin, 'omitnan'));
   end
   legend(currData.roiName, 'location', 'best', 'autoupdate', 'off')
   yL = ylim();
   for iTrial = 1:numel(currData.trialBounds{1})
       startVol = currData.trialBounds{1}(iTrial);
       plot([startVol, startVol], yL, '--', 'color', 'm')
   end
   
   
end

%% CHECK MEDIAN VALUES FOR PRE-ROI EXPERIMENTS

parentDir = 'F:\ImagingData';
expList = load_expList('groupName', 'singleTrialAcq_preROIs');

for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID)
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    
    blockFiles = dir(fullfile(parentDir, expDirName, 'sid_0', 'imagingData_reg_block*.mat'));
    medVals = [];
    for iFile = 1:numel(blockFiles)
        disp(iFile)
        load(fullfile(blockFiles(iFile).folder, blockFiles(iFile).name), 'imgData');
        sz = size(imgData);
        nVolumes = sz(4);
        rsData = reshape(imgData, sz(1), sz(2), sz(3), []);
        sz = size(rsData);
        medVals = [medVals, median(reshape(permute(rsData, [4 1 2 3]), sz(4), []), 2)'];        
    end
    expList.medVals{iExp} = medVals;
    expList.nVolumes(iExp) = nVolumes;
end

%% CHECK MEDIAN VALUES FOR PRE-FICTRAC 2018 EXPERIMENTS

expList = load_expList('groupName', 'singleTrialAcq_preFicTrac');
expList = expList(17:end, :)
    
    
for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID)
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    sidDir = dir(fullfile(parentDir, expDirName, 'sid_*'));
    if numel(sidDir) > 1
        blockFiles = dir(fullfile(parentDir, expDirName, 'sid_master', 'imagingData_reg_block*.mat'));
    else
        blockFiles = dir(fullfile(parentDir, expDirName, sidDir(1).name, 'imagingData_reg_block*.mat'));
    end
    
    medVals = [];
    for iFile = 1:numel(blockFiles)
        disp(iFile)
        try
        load(fullfile(blockFiles(iFile).folder, blockFiles(iFile).name), 'imgData');
        catch
            
        end
        sz = size(imgData);
        nVolumes = sz(4);
        rsData = reshape(imgData, sz(1), sz(2), sz(3), []);
        sz = size(rsData);
        medVals = [medVals, median(reshape(permute(rsData, [4 1 2 3]), sz(4), []), 2)'];        
    end
    expList.medVals{iExp} = medVals;
    expList.nVolumes(iExp) = nVolumes;
end
%% CHECK MEDIAN VALUES FOR PRE-FICTRAC 2017 EXPERIMENTS

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017';
expList = load_expList('groupName', 'singleTrialAcq_preFicTrac');
expList = expList(3:16, :)
    
    
for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID)
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    sidDir = dir(fullfile(parentDir, expDirName, 'sid_*'));
    if numel(sidDir) > 1
        blockFiles = dir(fullfile(parentDir, expDirName, 'sid_master', 'imagingData_reg_block*.mat'));
    else
        blockFiles = dir(fullfile(parentDir, expDirName, sidDir(1).name, 'imagingData_reg_block*.mat'));
    end
    
    medVals = [];
    for iFile = 1:numel(blockFiles)
        disp(iFile)
        try
        load(fullfile(blockFiles(iFile).folder, blockFiles(iFile).name), 'imgData');
        catch
            
        end
        sz = size(imgData);
        nVolumes = sz(4);
        rsData = reshape(imgData, sz(1), sz(2), sz(3), []);
        sz = size(rsData);
        medVals = [medVals, median(reshape(permute(rsData, [4 1 2 3]), sz(4), []), 2)'];        
    end
    expList.medVals{iExp} = medVals;
    expList.nVolumes(iExp) = nVolumes;
end
%%
for iExp = 1:size(expList, 1)
    
   % Create figure
   f = figure(iExp); clf; hold on
   f.Color = [1 1 1];
   f.Position = [50, 250, 1500, 700];
   plot(smoothdata(expList.medVals{iExp}, 'gaussian', 5));
   yL = ylim();
   trialBounds = 1:expList.nVolumes(iExp):numel(expList.medVals{iExp});
   for iTrial = 1:numel(trialBounds)
       startVol = trialBounds(iTrial);
       plot([startVol, startVol], yL, '--', 'color', 'm')
   end
end



















