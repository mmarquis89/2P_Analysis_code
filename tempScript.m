
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








