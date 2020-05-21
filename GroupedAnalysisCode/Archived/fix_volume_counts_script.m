% Fix uneven volume counts in gapless acquisition experiments


% Load exp list
expList = load_expList('groupName', 'gaplessAcq');


for iExp = 32:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID);
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    parentDir = find_parent_dir(currExpID);
    
    % Get list of unique sid/block combos
    dataFiles = dir(fullfile(parentDir, expDirName, 'cdata*.mat'));    
    fileNames = {dataFiles.name}';
    sidNums = str2double(regexp(fileNames, '(?<=sid_).', 'match', 'once'));
    bidNums = str2double(regexp(fileNames, '(?<=bid_).', 'match', 'once'));
    trialNums = str2double(regexp(fileNames, '(?<=_).....(?=\.mat)', 'match', 'once'));
    load(fullfile(parentDir, expDirName, 'volCounts.mat'), 'volCounts');
    trialList = table(sidNums, bidNums, trialNums, volCounts', fileNames, 'VariableNames', ...
            {'sid', 'blockNum', 'trialNum', 'volCounts', 'fileName'});
    blockList = unique(trialList(:, {'sid', 'blockNum'}));
    
    % Process each block individually
    for iBlock = 1:size(blockList, 1)
        
        currBlockTrials = innerjoin(blockList(iBlock, :), trialList);
        disp(['Sid ', num2str(currBlockTrials.sid(1)), ' block ' ...
            num2str(currBlockTrials.blockNum(1))])
        currVolCounts = currBlockTrials.volCounts;
        
        % Check whether total number of volumes fluctuates
        modeVols = mode(currVolCounts);
        remVols = unique(currVolCounts(currVolCounts ~= modeVols));
        if ~isnan(remVols)
            if numel(remVols) == 1 && abs(modeVols - remVols) < 2 
                
                % Trim one volume off all of the longer trials and re-save .mat file
                minVols = min([modeVols, remVols]);
                trimmedList = {};
                for iTrial = 1:size(currBlockTrials, 1)
                    currFileName = currBlockTrials.fileName{iTrial};
                    disp(currFileName);
                    load(fullfile(parentDir, expDirName, currFileName), 'imgData');
                    if size(imgData, 4) > minVols
                        trimmedList{end + 1} = currFileName;
                        imgData = imgData(:, :, :, 1:minVols, :);
                        save(fullfile(parentDir, expDirName, currFileName), 'imgData');
                    end
                end
                
                % Save a record of the files that were trimmed
                save(fullfile(parentDir, expDirName, ['sid_', num2str(currBlockTrials.sid(1)), ...
                        '_bid_', num2str(currBlockTrials.blockNum(1)), '_trimmedList.mat']), ...
                        'trimmedList');  
            else
                error('Volume count discrepancy')
            end
        end
    end%iBlock
end%iExp

