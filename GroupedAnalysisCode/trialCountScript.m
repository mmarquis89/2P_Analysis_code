
expList = load_expList('groupName', 'gaplessAcq');
imgDataNumbering = [];
uniqueBlockInfo = [];
for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID);
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    parentDir = find_parent_dir(currExpID);
    
    % Get a list of all the raw imaging data files in each gapless acquisition experiment
    dataFiles = dir(fullfile(parentDir, expDirName, 'cdata*.tif'));    
    fileNames = {dataFiles.name}';
    
    % Extract their sid, block numbers, and trial numbers
    sidNums = str2double(regexp(fileNames, '(?<=sid_).', 'match', 'once'));
    bidNums = str2double(regexp(fileNames, '(?<=bid_).', 'match', 'once'));
    trialNums = str2double(regexp(fileNames, '(?<=_).....(?=\.tif)', 'match', 'once'));
    
    % Save them all to a table with a row for each trial
    newRows = table(repmat({currExpID}, numel(sidNums), 1), sidNums, bidNums, trialNums, ...
            'VariableNames', {'expID', 'sid', 'blockNum', 'trialNum'});
    imgDataNumbering = [imgDataNumbering; newRows];
    
    % Also make a list of all unique sid and block numbers for each experiment
    uniqueSids = unique(sidNums);
    uniqueBids = unique(bidNums);
    newRow = table({currExpID}, {uniqueSids}, {uniqueBids}, 'VariableNames', {'expID', 'sids', ...
            'blockNums'});
    uniqueBlockInfo = [uniqueBlockInfo; newRow];  
   
end

% And finally one that gives the total number of trials in each block
uniqueBlocks = unique(imgDataNumbering(:, 1:3));
trialCounts = [];
maxTrialNums = [];
for iBlock = 1:size(uniqueBlocks, 1)
    currBlockTrials = innerjoin(uniqueBlocks(iBlock, :), imgDataNumbering);
    trialCounts(iBlock) = size(currBlockTrials, 1);
    maxTrialNums(iBlock) = currBlockTrials{end, end};
end
trialCountInfo = [uniqueBlocks, table(trialCounts', maxTrialNums', 'VariableNames', ...
        {'trialCount', 'maxTrialNum'})];
