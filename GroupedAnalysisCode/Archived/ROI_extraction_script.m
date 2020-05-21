% Load exp list
expList = load_expList('groupName', 'singleTrialAcq_preFicTrac');

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
    
    % Load summaryStats file (if one exists) for medVals offsetting
    if exist(fullfile(parentDir, expDirName, 'summaryStats.mat'), 'file')
        load(fullfile(parentDir, expDirName, 'summaryStats.mat'), 'summaryStats');
    end
    
    
    for iBlock = 1:numel(roiDefFiles)
        
        
    end
    
end%iExp