
saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';

% Load exp list
expList = load_expList('groupName', 'singleTrialAcq_preROIs');
expList = [expList; load_expList('groupName', 'singleTrialAcq_with_ROIs')];

for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID);
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    parentDir = find_parent_dir(currExpID);
    
    % Skip experiment if there are any ROI def files in the directory
    roiDefFiles = dir(fullfile(parentDir, expDirName, 'roiDefs_*.mat')); 
    if isempty(roiDefFiles)
        disp(['No ROI def files found for ', expDirName, '...skipping experiment']);
        continue 
    end
    nBlocks = numel(roiDefFiles);
    
    % Load summaryStats file (if one exists) for medVals offsetting
    if exist(fullfile(parentDir, expDirName, 'summaryStats.mat'), 'file')
        medValOffset = 1;
        load(fullfile(parentDir, expDirName, 'summaryStats.mat'), 'medVals');
    else
        medValOffset = 0;
    end
    
    trialCount = 0;
    roiData = [];
    blockMetadata = [];
    for iBlock = 1:nBlocks
        disp(['Block ', num2str(iBlock), ' of ', num2str(nBlocks)])
        
        % Load imgData and roiDefs for current block
        disp('Loading data...')
        load(fullfile(parentDir, expDirName, ['imagingData_reg_block_', ...
                pad(num2str(iBlock), 3, 'left', '0'), '.mat']), 'imgData');
        load(fullfile(parentDir, expDirName, ['roiDefs_block_', ...
                pad(num2str(iBlock), 3, 'left', '0'), '.mat']), 'roiDefs');
        
        % Process each trial individually
        nBlockTrials = size(imgData, 5);
        if medValOffset
            currMedVals = medVals([1:nBlockTrials] + trialCount);
        end
        disp('Extracting data from ROIs...');
        clipTotals = zeros(numel(roiDefs), 1);
        for iTrial = 1:nBlockTrials
            
            % Subtract median offset, if available
            if medValOffset
                trialMedVal = median(currMedVals{iTrial});
            else
                trialMedVal = 0;
            end
            currTrialData = imgData(:, :, :, :, iTrial) - trialMedVal; % --> [y, x, plane, volume]
            
            % Reshape into 1D frames
            sz = size(currTrialData);
            rsData = reshape(permute(currTrialData, [3 4 1 2]), sz(3), sz(4), []); % --> [plane, volume, pixel]
            
            % Loop through and extract data for each subROI
            for iROI = 1:numel(roiDefs)
                currROI = roiDefs(iROI);
                
                % Save block metadata on first trial of each block
                if iTrial == 1
                    newRow = table({currExpID}, iBlock, {currROI.name}, {currROI.subROIs}, ...
                            {(1:nBlockTrials) + trialCount},...
                            'VariableNames', {'expID', 'blockNum', 'roiName', 'subROIs', ...
                            'trialNums'});
                    if medValOffset
                        newRow.medVals = {currMedVals};
                    else
                        newRow.medVals = {{}};
                    end
                    blockMetadata = [blockMetadata; newRow];
                end
                
                currRoiData = [];
                for iSubROI = 1:numel(currROI.subROIs)
                    currSubROI = currROI.subROIs(iSubROI);
                    currImgData = squeeze(rsData(currSubROI.plane, :, :));  % --> [volume, pixel]
                    mask = currSubROI.mask(:);                              % --> [pixel]
                    currImgData = currImgData(:, mask)';                    % --> [pixel, volume]
                    currRoiData = [currRoiData; currImgData];               % --> [pixel, volume]
                end
                
                % Average data across pixels and subROIs
                roiDataAvg = mean(currRoiData, 1); % --> [volume]
                
                % Create new row and append to ROI data table
                newRow = table({currExpID}, 'variableNames', {'expID'});
                newRow.trialNum = iTrial + trialCount;
                newRow.roiName = {roiDefs(iROI).name};
                newRow.subROIs = {rmfield(roiDefs(iROI).subROIs, {'obj', 'refImage', 'mask'})};
                newRow.rawFl = {roiDataAvg'};
                roiData = [roiData; newRow];
                
            end%iROI
        end%iTrial
        
        % Increment trial count
        trialCount = trialCount + nBlockTrials;
        
    end%iBlock
    disp('All ROI data extracted')
    
    % Subtract minimum value across entire experiment from all ROI data if it is < 0
    expMinVal = min(cell2mat(roiData.rawFl));
    if expMinVal < 0
        disp(['Subtracting minimum value of ', num2str(expMinVal) ' from all ROI data'])
        for iRow = 1:size(roiData, 1)
            roiData.rawFl{iRow} = (roiData.rawFl{iRow} - expMinVal) + 1; % Add one to allow division
            
        end
    end
    
    % Calculate trial-based dF/F, using median of bottom 5th percentile as baseline
    disp('Calculating dF/F...')
    dffData = [];
    for iRow = 1:size(roiData, 1)
        roiDataSorted = sort(roiData.rawFl{iRow});
        baselineF = median(roiDataSorted(1:round(numel(roiDataSorted) * 0.05)));
        dffData{iRow} = (roiData.rawFl{iRow} - baselineF) ./ baselineF;
    end
    roiData = [roiData, table(dffData', 'variableNames', {'dffData'})];
    
    % Calculate dF/F using median of bottom 5th percentile across entire experiment as a baseline
    roiList = unique(roiData(:, 3));
    baselineVals = [];
    for iROI = 1:size(roiList, 1)
        currRoiData = innerjoin(roiList(iROI, :), roiData);
        roiDataSort = sort(cell2mat(currRoiData.rawFl));
        baselineVals(iROI) = median(roiDataSort(1:round(numel(roiDataSort) * 0.05)));        
    end
    baselineTable = [roiList, table(baselineVals', 'variableNames', {'baselineVals'})];
    expDffData = [];
    for iRow = 1:size(roiData, 1)
        currRowJoin = innerjoin(roiData(iRow, :), baselineTable);
        expDffData{iRow} = (currRowJoin.rawFl{:} - currRowJoin.baselineVals) / currRowJoin.baselineVals;
    end
    roiData = [roiData, table(expDffData', 'variableNames', {'expDffData'})];
    
    % Save current experiment's ROI data
    disp('Saving ROI data...')
    save(fullfile(saveDir, [currExpID, '_roiData.mat']), 'roiData')
     
    % Save subROI data separately to save space
    save(fullfile(saveDir, [currExpID, '_blockMetadata.mat']), 'blockMetadata');
    
    
    
end%iExp


















