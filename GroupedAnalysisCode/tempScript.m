
% Create DataTable with metadata and ROI data
expList = load_expList;
expList = expList(contains(expList.expName, 'D-ANT'), :);
expList = expList(~contains(expList.groupName, 'newExpts'), :);
expList(20, :) = [];
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';

allExpData = [];

for iExp = 1:size(expList, 1)
    
    currExpID = expList.expID{iExp};
    disp(currExpID);
    
    % Generate full paths to data files
    roiFileName = [currExpID, '_roiData.mat'];
    roiFile = fullfile(parentDir, roiFileName);
    ftFileName = [currExpID, '_ficTracData.mat'];
    ftFile = fullfile(parentDir, ftFileName);
    expMdFile = fullfile(parentDir, [currExpID, '_expMetadata.csv']);
    trialMdFile = fullfile(parentDir, [currExpID, '_trialMetadata.mat']);
    
    if exist(roiFile, 'file') && exist(ftFile, 'file') && exist(expMdFile, 'file')
        load(roiFile, 'roiData');
        load(ftFile, 'ftData');
        expMd = readtable(expMdFile, 'delimiter', ',');
        load(trialMdFile, 'trialMetadata');
        trialMd = trialMetadata;
        if contains(expList.groupName{iExp}, 'gapless')
            acqNumTable = trialMd(:, {'expID', 'trialNum'});
            acqNumTable.acqNum = acqNumTable.trialNum;
        else 
            load(fullfile(parentDir, [currExpID, '_blockMetadata.mat']), 'blockMetadata');            
            blockMetadata = blockMetadata(strcmpi(blockMetadata.roiName, 'Background'), :);
            trialCounts = blockMetadata.trialNums;
            acqNumTable = trialMd(:, {'expID', 'trialNum'});
            acqNums = ones(size(acqNumTable, 1), 1);
            for iAcq = 1:numel(trialCounts) 
                acqNums(trialCounts{iAcq}, 1) = blockMetadata.blockNum(iAcq);
            end
            acqNumTable.acqNum = acqNums;
        end
        allExpData = [allExpData; inner_join(expMd, acqNumTable, trialMd, roiData, ftData)];

    else
        disp(['Skipping ', currExpID, ' due to one or more missing files']);
    end
end%iExp

% Create source DataTable
sourceDataTable = DataTable(allExpData);

%%
expNum = 27;

filterDefs = struct();
filterDefs.roiName = 'TypeD';
plot_1_maxFl = 600;
% filterDefs.roiName = 'ANT';
% plot_1_maxFl = 800;
plot_1_maxSpeed = 20;
plot_2_max = 150;
smWinVols = 3;
smWinFrames = 3;
smReps = 10;
spacerRows = 2;

filterDefs.expID = expList.expID{expNum};
disp(filterDefs.expID)

sourceDataTable = sourceDataTable.initialize_filters(filterDefs);
uniqueROIs = unique(sourceDataTable.subset.roiName);
if numel(uniqueROIs) > 1
    filterDefs.roiName = uniqueROIs{1};
    sourceDataTable = sourceDataTable.initialize_filters(filterDefs);
    disp('Bilateral ROIs detected...using first ROI only')
end
subsetData = sourceDataTable.subset;
disp(unique(subsetData.roiName))

FRAME_RATE = 25;
blockFlMats = []; blockSpeedMats = []; stimOnsets = []; stimOffsets = [];
if all(subsetData.trialNum == subsetData.acqNum)
    
    % ----- Gapless acqusition -----
    for iBlock = 1:numel(unique(subsetData.trialNum))
        currBlockData = subsetData(subsetData.acqNum == iBlock, :);
        flData = currBlockData.rawFl{:};
        if ~isempty(currBlockData.pmtShutoffVols{:}) & ~isnan(currBlockData.pmtShutoffVols{:}) %#ok<AND2>
            flData(currBlockData.pmtShutoffVols{:}) = nan;
        end
        blockFlMats{iBlock} = reshape(flData, [], currBlockData.originalTrialCount);
        blockSpeedMats{iBlock} = reshape(currBlockData.moveSpeed{:}, [], ...
                currBlockData.originalTrialCount);
        
        % Find stim onsets and offsets
        volTimes = calc_volTimes(numel(flData), currBlockData.volumeRate, ...
                currBlockData.trialDuration, currBlockData.originalTrialCount);
        allEvents = load_event_data(expList.expID(expNum));
        eventNames = fieldnames(allEvents);
        allStimVec = zeros(size(volTimes));
        for iName = 1:numel(eventNames)
            currName = eventNames{iName};
            if ismember(currName, {'odor', 'soundstim', 'optostim'})                
                eventArr = allEvents.(currName).create_logical_array(numel(volTimes), volTimes, ...
                    currBlockData(:, {'expID', 'trialNum'}));
                allStimVec = allStimVec | eventArr';
            end
        end
        currStimOnsets = zeros(size(volTimes));
        currStimOffsets = zeros(size(volTimes));
        currStimOnsets(strfind(allStimVec, [0 1])) = 1;
        currStimOffsets(strfind(allStimVec, [1 0])) = 1;
        if allStimVec(end) == 1
            currStimOffsets(end) = 1;
        end
        if allStimVec(1) == 1;
           currStimOnsets(1) = 1; 
        end
        stimOnsets{iBlock} = reshape(currStimOnsets, [], currBlockData.originalTrialCount);
        stimOffsets{iBlock} = reshape(currStimOffsets, [], currBlockData.originalTrialCount);
    end    
else
    % ----- Single-trial acquistion -----
    flMat = cell2padded_mat(subsetData.rawFl); % --> [vol, trial]
    moveSpeedMat = cell2padded_mat(subsetData.moveSpeed); % --> [frame, trial]
        
    % Set PMT shutoff vols to NaN
    for iTrial = 1:size(flMat, 2)
        if ~isempty(subsetData.pmtShutoffVols{iTrial}) & ~isnan(subsetData.pmtShutoffVols{iTrial}) %#ok<AND2>
            flMat(subsetData.pmtShutoffVols{iTrial}, iTrial) = nan;
        end
    end
    
    % Find stim onsets and offsets
    volTimes = calc_volTimes(size(flMat, 1), subsetData.volumeRate(1), ...
            subsetData.trialDuration(1), 1);
    allEvents = load_event_data(expList.expID(expNum));
    eventNames = fieldnames(allEvents);
    allStimArr = zeros(size(flMat));
    for iName = 1:numel(eventNames)
        currName = eventNames{iName};
        if ismember(currName, {'odor', 'soundstim', 'optostim'}) 
            currEventArr = allEvents.(currName).create_logical_array(numel(volTimes), volTimes, ...
                subsetData(:, {'expID', 'trialNum'}));
            allStimArr = allStimArr | currEventArr;
        end
    end
    allStimOnsets = zeros(size(flMat));
    allStimOffsets = zeros(size(flMat));
    for iTrial = 1:size(allStimArr, 2)
        allStimOnsets(strfind(allStimArr(:, iTrial)', [0 1]), iTrial) = 1;
        allStimOffsets(strfind(allStimArr(:, iTrial)', [1 0]), iTrial) = 1;
    end
    
    % Separate into one block for each acqusition
    for iBlock = 1:numel(unique(subsetData.acqNum))
        blockFlMats{iBlock} = flMat(:, subsetData.acqNum == iBlock);
        blockSpeedMats{iBlock} = moveSpeedMat(:, subsetData.acqNum == iBlock);        
        stimOnsets{iBlock} = allStimOnsets(:, subsetData.acqNum == iBlock);
        stimOffsets{iBlock} = allStimOffsets(:, subsetData.acqNum == iBlock);
    end
end

f = figure(1);clf; 
f.Color = [1 1 1];
    %     if strcmp(flType, 'trialDff')
    %         baseF = repmat(obj.sourceDataTable.subset.trialBaseline', size(outputMat, 1), 1);
    %         outputMat = (outputMat - baseF) ./ baseF;
    %     elseif strcmp(flType, 'expDff')
    %         baseF = repmat(obj.sourceDataTable.subset.expBaseline', size(outputMat, 1), 1);
%         outputMat = (outputMat - baseF) ./ baseF;
%     end
flMat = cell2mat(blockFlMats)';                             % --> [trial, volume]
speedMat = cell2mat(blockSpeedMats)' .* 4.5 .* FRAME_RATE;   % --> [trial, frame]
smFlMat = smoothdata(flMat, 2, 'gaussian', smWinVols);
smSpeedMat = repeat_smooth(speedMat, smReps, 'dim', 2, 'smWin', smWinFrames);
smFlMat(smFlMat > plot_1_maxFl) = plot_1_maxFl;
smSpeedMat(smSpeedMat > plot_1_maxSpeed) = plot_1_maxSpeed;

trialGroups = repelem(1:numel(blockFlMats), cellfun(@(x) size(x, 2), blockFlMats));
ax = subaxis(1, 2, 1);
plot_2D_summary(smFlMat, subsetData.volumeRate(1), 'plotaxes', ax, 'trialGroups', trialGroups, ...
        'spacerArrRows', spacerRows);
colorbar()
hold on;
yCount = 0.5;
yL = ylim();
for iBlock = 1:numel(stimOnsets)
    currOnsets = stimOnsets{iBlock};
    currOffsets = stimOffsets{iBlock};
    for iTrial = 1:size(currOnsets, 2)
        currOnsetVols = find(currOnsets(:, iTrial));
        currOffsetVols = find(currOffsets(:, iTrial));
        yy = [yCount, yCount + 1];
        for iStim = 1:numel(currOnsetVols)
            xx = [1 1] * currOnsetVols(iStim);
            plot(xx, yy, 'linewidth', 2, 'color', 'g');
            xx = [1 1] * currOffsetVols(iStim);
            plot(xx, yy, 'linewidth', 2, 'color', 'r');
        end
        yCount = yCount + 1;
    end
    yCount = yCount + spacerRows;
end
ylim(yL);

title('Raw F');

ax = subaxis(1, 2, 2);
plot_2D_summary(smSpeedMat, FRAME_RATE, 'plotaxes', ax, 'trialGroups', trialGroups, ...
        'spacerArrRows', 2);
colorbar();
hold on;
vol2frame = sample_lookup(FRAME_RATE, subsetData.volumeRate(1));
yCount = 0.5;
yL = ylim();
for iBlock = 1:numel(stimOnsets)
    currOnsets = stimOnsets{iBlock};
    currOffsets = stimOffsets{iBlock};
    for iTrial = 1:size(currOnsets, 2)
        if any(currOnsets(:, iTrial))
            currOnsetFrames = vol2frame.convert(find(currOnsets(:, iTrial)));
            currOffsetFrames = vol2frame.convert(find(currOffsets(:, iTrial)));
            yy = [yCount, yCount + 1];
            for iStim = 1:numel(currOffsetFrames)
                xx = [1 1] * currOnsetFrames(iStim);
                plot(xx, yy, 'linewidth', 2, 'color', 'g');
                xx = [1 1] * currOffsetFrames(iStim);
                plot(xx, yy, 'linewidth', 2, 'color', 'r');
            end
        end
        yCount = yCount + 1;
    end
    yCount = yCount + spacerRows;
end
ylim(yL);
title('Move speed (mm/sec)')

%%
filterDefs.roiName = 'TypeB';
sourceDataTable = sourceDataTable.initialize_filters(filterDefs);
subsetData = sourceDataTable.subset;

disp(filterDefs.expID)
disp(unique(sourceDataTable.subset.roiName))

nChunksPerBlock = [];
chunkedFl = [];
if size(subsetData, 1) < 20
    
    % Gapless acquisition
    for iTrial = 1:size(subsetData, 1)
        currFl = subsetData.rawFl{iTrial};
        if ~isempty(subsetData.pmtShutoffVols{iTrial}) && ~any(isnan(subsetData.pmtShutoffVols{iTrial}))
            currFl(subsetData.pmtShutoffVols{iTrial}) = nan;
        end
        nChunks = round((numel(currFl) / subsetData.volumeRate(iTrial)) / chunkLength);
        chunkSize = floor(numel(currFl) / nChunks);
        extraVols = numel(currFl) - (nChunks * chunkSize);
        while extraVols > 0
            nChunks = nChunks + 1;
            nanVols = (nChunks * chunkSize) - numel(currFl);
            currFl = [currFl; nan(nanVols, 1)];
            extraVols = numel(currFl) - (nChunks * chunkSize);
        end
        rsFl = reshape(currFl, [], nChunks);
        chunkedFl = [chunkedFl, rsFl];
        nChunksPerBlock(iTrial) = size(rsFl, 2);
    end
else
    % Single-trial acquisition
    currSubset = sourceDataTable.subset;
    chunkedFl = [chunkedFl, cell2padded_mat(currSubset.rawFl)]; % --> [vol, trial]
    for iTrial = 1:size(chunkedFl, 2)
        if ~isempty(subsetData.pmtShutoffVols{iTrial})
            chunkedFl(subsetData.pmtShutoffVols{iTrial}, iTrial) = nan;
        end
    end
    nChunksPerBlock = size(chunkedFl, 2);
    %     if strcmp(flType, 'trialDff')
    %         baseF = repmat(obj.sourceDataTable.subset.trialBaseline', size(outputMat, 1), 1);
    %         outputMat = (outputMat - baseF) ./ baseF;
    %     elseif strcmp(flType, 'expDff')
    %         baseF = repmat(obj.sourceDataTable.subset.expBaseline', size(outputMat, 1), 1);
%         outputMat = (outputMat - baseF) ./ baseF;
%     end
    
end
chunkedFl(chunkedFl > plot_2_max) = plot_2_max;
subaxis(1, 2, 2);
imagesc(chunkedFl')
colorbar();
title('TypeB')

















