
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\gaplessAcq';

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
allExpOdorEvents = odorEvent();
roiDataTable = [];
ftDataTable = [];
for iExp = 1:numel(gaplessAcqExpList)
    
    currExpID = gaplessAcqExpList{iExp};
    disp(currExpID);
    
    % trialMetadata files
    currTrialMd = readtable(fullfile(parentDir, [currExpID, '_trialMetadata.csv']));
    allTrialMd = [allTrialMd; currTrialMd];
    
    % Odor event data
    odorEventFile = fullfile(parentDir, [currExpID, '_event_data_odor.csv']);
    if exist(odorEventFile, 'file')
        allExpOdorEvents = allExpOdorEvents.load_csv(parentDir, 'fileNamePrefix', currExpID);
    end
    
    % ROI data
    load(fullfile(parentDir, [currExpID, '_roiData.mat']), 'roiData');
    for iTrial = 1:numel(roiData)
        for iROI = 1:numel(roiData{iTrial})
            newRow = table({currExpID}, 'VariableNames', {'expID'});
            newRow.trialNum = iTrial;
            newRow.roiName = {roiData{iTrial}(iROI).name};
            newRow.rawFl = {roiData{iTrial}(iROI).rawFl};
            newRow.expDffData = {roiData{iTrial}(iROI).expDffData};
            roiDataTable = [roiDataTable; newRow];
        end
    end% iTrial
    
%     % FicTrac data
%     load(fullfile(parentDir, [currExpID, '_ficTracData.mat']), 'ftData');    
%     newData = [table(repmat({currExpID}, numel(ftData), 1), 'VariableNames', {'expID'}), ...
%             struct2table(ftData, 'AsArray', 1)];
%     ftDataTable = [ftDataTable; newData];
%     
end%iExp


%% Make list of number of experiments for each unique odor identity/concentration combo

uniqueStims = unique(allExpOdorEvents(:, 5:6));
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
analysisWin = [8 12];

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
            stimCount = stimCount + 1;
        end
    end
    outputTable.expDff{iExp} = currExpDff;
end

%% Fix ROI data


saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\gaplessAcq';
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData'; 
expList = readtable(fullfile(parentDir, 'oldExpList_gaplessAcq.csv'), 'delimiter', ',');
expList = expList([1:10, 12:end],:); % Because I haven't done 11/9 #1 yet

for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID)
    
    % Load trial metadata
    trialMd = readtable(fullfile(saveDir, [currExpID, '_trialMetadata.csv']), 'delimiter', ',');
    
    % Load ROI data
    currFile = fullfile(saveDir, [currExpID, '_roiData.mat']);
    load(currFile, 'roiData');
    figure(iExp);clf;hold on
    for iROI = 1:size(roiData, 1)
        
        plot(roiData{iROI, 'rawFl'}{:})
        title(currExpID)
        
        
%         currRoiMd = innerjoin(roiData(iROI, :), trialMd);
%         trialCount = currRoiMd.originalTrialCount;
%         nVolumes = currRoiMd.nVolumes / trialCount;            
%             
%         roiData{iROI, 'rawFl'}{:} = as_vector(reshape(roiData{iROI, 'rawFl'}{:}, trialCount, ...
%                 nVolumes)');
%         roiData{iROI, 'dffData'}{:} = as_vector(reshape(roiData{iROI, 'dffData'}{:}, trialCount, ...
%                 nVolumes)');
%         roiData{iROI, 'expDffData'}{:} = as_vector(reshape(roiData{iROI, 'expDffData'}{:}, ...
%                 trialCount, nVolumes)');
    end
    
    % Save updated ROI data
%     save(currFile, 'roiData');
    
end









