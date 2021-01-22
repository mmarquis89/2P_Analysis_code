
saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
% 
% % Load exp list
% expList = load_expList('groupName', 'gaplessAcq');

expList = table({'20190315-3'}, 'variablenames', {'expID'}); 
expList.expName = {'D-ANT_6s'};

for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(' ')
    disp('**********')
    disp(currExpID);
    disp('**********')
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    parentDir = find_parent_dir(currExpID);
    
    % Skip experiment if there aren't any ROI def files in the directory
    roiDefFiles = dir(fullfile(parentDir, expDirName, 'roiDefs_*.mat')); 
    if isempty(roiDefFiles)
        disp(['No ROI def files found for ', expDirName, '...skipping experiment']);
        continue 
    end
    
    acqNums = cellfun(@str2double, unique(regexp({roiDefFiles.name}, '(?<=block_)...(?=-.\.mat)', ...
            'match', 'once')));
    nAcqs = numel(acqNums);
    roiData = [];
    for iAcq = 1:nAcqs
        blockFiles = dir(fullfile(parentDir, expDirName, ['roiDefs_block_', pad(num2str(iAcq), ...
                3, 'left', '0'), '*.mat']));

        for iBlock = 1:numel(blockFiles)
            disp(['Acq #', num2str(iAcq), ' block #', num2str(iBlock)])
            
            % Load imgData and roiDefs for current block
            disp('Loading data...')
            load(fullfile(parentDir, expDirName, regexprep(blockFiles(iBlock).name, ...
                    'roiDefs_', 'imagingData_reg_')), 'imgData');
            load(fullfile(parentDir, expDirName, blockFiles(iBlock).name), 'roiDefs');
            
            % Reshape imaging data to get rid of artifical trial boundaries
            sz = size(imgData);
            rsImgData = reshape(imgData, sz(1), sz(2), sz(3), []); % --> [y, x, plane, volume]
            
            % Then reshape again into 1D frames
            sz = size(rsImgData);
            rsData = reshape(permute(rsImgData, [3 4 1 2]), sz(3), sz(4), []); % --> [plane, volume, pixel]
            
            % Create table to hold data for current acquisition
            if iBlock == 1
                currAcqRoiData = [];
                for iROI = 1:numel(roiDefs)
                    newRow = table({currExpID}, iAcq, {roiDefs(iROI).name}, ...
                        {roiDefs(iROI).subROIs}, {[]}, ...
                        'VariableNames', {'expID', 'trialNum', 'roiName', 'subROIs', 'rawFl'});
                    currAcqRoiData = [currAcqRoiData; newRow];
                end
                
            end
            
            % Loop through and extract data for each subROI
            for iROI = 1:numel(roiDefs)
                currROI = roiDefs(iROI);
                currRoiData = [];
                for iSubROI = 1:numel(currROI.subROIs)
                    currSubROI = currROI.subROIs(iSubROI);
                    currImgData = squeeze(rsData(currSubROI.plane, :, :));  % --> [volume, pixel]
                    mask = currSubROI.mask(:);                              % --> [pixel]
                    currImgData = currImgData(:, mask)';                    % --> [pixel, volume]
                    currRoiData = [currRoiData; currImgData];               % --> [pixel, volume]
                end
                
                % Average data across pixels and subROIs
                roiDataAvg = mean(currRoiData, 1)'; % --> [volume]  
                currAcqRoiData.rawFl{iROI} = [currAcqRoiData.rawFl{iROI}; roiDataAvg];
                
            end%iROI
        end%iBlock
        
        % Append current acquisition's data to the table after concatenating all blocks
        roiData = [roiData; currAcqRoiData];
        
    end%iAcq
    
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
%      
%     % Save subROI data separately to save space
%     save(fullfile(saveDir, [currExpID, '_blockMetadata.mat']), 'blockMetadata');
    
    
    
end%iExp


















