
% Create DataTable with metadata and ROI data
expList = load_expList;
expList = expList(contains(expList.expName, 'PPM2'), :);
expList = expList(~contains(expList.groupName, 'newExpts'), :);
% expList(20, :) = [];
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
expNum = 11;

roiNames = {'TypeF', 'TypeF'};
% roiNames = {'TypeD', 'TypeB'};

ROI_1_maxFl = 600;
ROI_2_maxFl = 600;
maxSpeed = 25;

smWinVols = 3;
smWinFrames = 3;
smReps = 10;
spacerRows = 2;

SV = 0.05;
SH = 0.03;
ML = 0.04;
MR = 0.015;
MT = 0.04;
MB = 0.09;

filterDefs = struct();
filterDefs.expID = expList.expID{expNum};
disp(filterDefs.expID)

% Extract primary ROI data
filterDefs.roiName = roiNames{1};
sourceDataTable = sourceDataTable.initialize_filters(filterDefs);
uniqueROIs = unique(sourceDataTable.subset.roiName);
if numel(uniqueROIs) > 1
    filterDefs.roiName = uniqueROIs{1};
    sourceDataTable = sourceDataTable.initialize_filters(filterDefs);
    disp('Bilateral ROIs detected...using first ROI only')
end
subsetData = sourceDataTable.subset;
disp(unique(subsetData.roiName))

% Extract comparison ROI data
filterDefs.roiName = roiNames{2};
sourceDataTable = sourceDataTable.initialize_filters(filterDefs);
uniqueROIs = unique(sourceDataTable.subset.roiName);
if numel(uniqueROIs) > 1
    filterDefs.roiName = uniqueROIs{1};
    sourceDataTable = sourceDataTable.initialize_filters(filterDefs);
    disp('Bilateral ROIs detected...using first ROI only')
end
subsetData_compROI = sourceDataTable.subset;
disp(unique(subsetData_compROI.roiName))

FRAME_RATE = 25;
blockFlMats = []; blockFlMats_compROI = []; blockSpeedMats = []; stimOnsets = []; stimOffsets = [];
if all(subsetData.trialNum == subsetData.acqNum)
    
    % ----- Gapless acqusition -----
    for iBlock = 1:numel(unique(subsetData.trialNum))
        currBlockData = subsetData(subsetData.acqNum == iBlock, :);
        currBlockData_compROI = subsetData_compROI(subsetData_compROI.acqNum == iBlock, :);
        flData = currBlockData.rawFl{:};
        flData_compROI = currBlockData_compROI.rawFl{:};
        if ~isempty(currBlockData.pmtShutoffVols{:}) & ~isnan(currBlockData.pmtShutoffVols{:})
            flData(currBlockData.pmtShutoffVols{:}) = nan;
            flData_compROI(currBlockData_compROI.pmtShutoffVols{:}) = nan;
        end
        blockFlMats{iBlock} = reshape(flData, [], currBlockData.originalTrialCount);
        blockFlMats_compROI{iBlock} = reshape(flData_compROI, [], ...
                currBlockData_compROI.originalTrialCount);
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
        if allStimVec(1) == 1
           currStimOnsets(1) = 1; 
        end
        stimOnsets{iBlock} = reshape(currStimOnsets, [], currBlockData.originalTrialCount);
        stimOffsets{iBlock} = reshape(currStimOffsets, [], currBlockData.originalTrialCount);
    end    
else
    % ----- Single-trial acquistion -----
    flMat = cell2padded_mat(subsetData.rawFl); % --> [vol, trial]
    flMat_compROI = cell2padded_mat(subsetData_compROI.rawFl); % --> [vol, trial]
    moveSpeedMat = cell2padded_mat(subsetData.moveSpeed); % --> [frame, trial]
        
    % Set PMT shutoff vols to NaN
    for iTrial = 1:size(flMat, 2)
        if ~isempty(subsetData.pmtShutoffVols{iTrial}) & ~isnan(subsetData.pmtShutoffVols{iTrial}) %#ok<AND2>
            flMat(subsetData.pmtShutoffVols{iTrial}, iTrial) = nan;
            flMat_compROI(subsetData_compROI.pmtShutoffVols{iTrial}, iTrial) = nan;
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
    
    % Treating individual trial acquisitions as one big block
    blockFlMats = {flMat};
    blockFlMats_compROI = {flMat_compROI};
    blockSpeedMats = {moveSpeedMat};
    stimOnsets = {allStimOnsets};
    stimOffsets = {allStimOffsets};
    
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

% Process data matrices
flMat = cell2mat(blockFlMats)';                             % --> [trial, volume]
flMat_compROI = cell2mat(blockFlMats_compROI)';                             % --> [trial, volume]
speedMat = cell2mat(blockSpeedMats)' .* 4.5 .* FRAME_RATE;   % --> [trial, frame]
smFlMat = smoothdata(flMat, 2, 'gaussian', smWinVols);
smFlMat_compROI = smoothdata(flMat_compROI, 2, 'gaussian', smWinVols);
smSpeedMat = repeat_smooth(speedMat, smReps, 'dim', 2, 'smWin', smWinFrames);
smFlMat(smFlMat > ROI_1_maxFl) = ROI_1_maxFl;
smFlMat_compROI(smFlMat_compROI > ROI_2_maxFl) = ROI_2_maxFl;
smSpeedMat(smSpeedMat > maxSpeed) = maxSpeed;

trialGroups = repelem(1:numel(blockFlMats), cellfun(@(x) size(x, 2), blockFlMats));

% Create fluorescence plot for primary (Type D) ROI
ax = subaxis(1, 3, 1, 'sv', SV, 'sh', SH, 'ml', ML, 'mr', MR, 'mt', MT, 'mb', MB);
plot_2D_summary(smFlMat, subsetData.volumeRate(1), 'plotaxes', ax, 'trialGroups', trialGroups, ...
        'spacerArrRows', spacerRows);
colorbar(); hold on;
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
title([roiNames{1}, '  -  Raw F']);

% Create moveSpeed plot
ax = subaxis(1, 3, 2);
plot_2D_summary(smSpeedMat, FRAME_RATE, 'plotaxes', ax, 'trialGroups', trialGroups, ...
        'spacerArrRows', 2);
colorbar(); hold on;
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
ax.YLabel.String = '';

% Create fluorescence plot for comparison ROI
ax = subaxis(1, 3, 3);
plot_2D_summary(smFlMat_compROI, subsetData.volumeRate(1), 'plotaxes', ax, 'trialGroups', trialGroups, ...
        'spacerArrRows', spacerRows);
colorbar(); hold on;
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
title([roiNames{2}, '  -  Raw F']);
ax.YLabel.String = '';
