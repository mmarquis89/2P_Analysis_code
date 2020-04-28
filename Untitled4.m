saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';

expList = readtable(fullfile(saveDir, 'oldExpList_preFicTrac.csv'), 'delimiter', ',');

allFilesList = [];
dirInfoTable = []; dirCheckList = [];
for iExp = 1:size(expList, 1)
   
    currExpID = expList.expID{iExp};
    disp(currExpID);
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    
    % Find experiment directory  
    if strcmp(expDirName(1:4), '2017')
        parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017'; 
    else
        parentDir = 'F:\ImagingData'; 
    end
    sidDir = dir(fullfile(parentDir, expDirName, 'sid_*'));
    sidDir = sidDir([sidDir.isdir]);
    if numel(sidDir) == 1
        expDir = fullfile(parentDir, expDirName, sidDir.name); % If there's only one sid directory
    elseif numel(sidDir) > 1
        expDir = fullfile(parentDir, expDirName, 'sid_master'); % If I've already manually combined sids
    else
        error('Could not find experiment directory');
    end
    expName = expList.expName{iExp};
    
    expDirFiles = dir(expDir);
    allFilesList = [allFilesList; table({expDirFiles.name}')];
    for iFile = 1:numel(expDirFiles)
        newRow =  table({currExpID}, {expDirFiles(iFile).name});
        dirInfoTable = [dirInfoTable; newRow];
    end
     
    file_exists = @(x) any(~cellfun(@isempty, regexp(dirInfoTable{strcmp(dirInfoTable.Var1, ...
            currExpID), 2}, x)));
    
    dirCheckList(iExp, 1) = file_exists('imagingData_reg_block_..\.mat');
    dirCheckList(iExp, 2) = file_exists('refImages_reg_block_..\.mat');
    dirCheckList(iExp, 3) = file_exists('refImages_sid_._sessionFile\.mat');
    dirCheckList(iExp, 4) = file_exists('sid_._refImages\.mat');
    dirCheckList(iExp, 5) = file_exists('rigid_sid_._sessionFile\.mat');
    dirCheckList(iExp, 6) = file_exists('rigid_sid_._sessionFile_Reg1\.mat');
    dirCheckList(iExp, 7) = file_exists('rigid_sid_._Chan_1_sessionFile\.mat');
    dirCheckList(iExp, 8) = file_exists('rigid_sid_._Chan_1_sessionFile_Reg1\.mat');
    dirCheckList(iExp, 9) = file_exists('sid_._Chan_1_sessionFile_Reg1\.mat');
    dirCheckList(iExp, 10) = file_exists('sid_._Annotations\.mat');
    dirCheckList(iExp, 11) = file_exists('sid_._BehavioralAnnotations\.mat');
    dirCheckList(iExp, 12) = file_exists('sid_._RawFrames\.7z');
    dirCheckList(iExp, 13) = file_exists('metadata_.*\.mat');
    dirCheckList(iExp, 14) = file_exists('bdata_.*\.mat');
    dirCheckList(iExp, 15) = file_exists('cdata_.*\.tif');
    dirCheckList(iExp, 16) = file_exists('imagingData_raw.*\.mat');
    dirCheckList(iExp, 17) = file_exists('imagingData_reg.*\.mat');

end%iExp

dirCheckList = logical(dirCheckList);
checkListTable = array2table(dirCheckList, 'VariableNames', {'imagingData_reg_block_X', ...
        'refImages_reg_block_X', 'refImages_sid_X_sessionFile', ...
        'sid_X_refImages', 'rigid_sid_X_sessionFile', 'rigid_sid_X_sessionFile_Reg1', ...
        'rigid_sid_X_Chan_1_sessionFile', 'rigid_sid_X_Chan_1_sessionFile_Reg1', ...
        'sid_X_Chan_1_sessionFile_Reg1', 'sid_X_Annotations', 'sid_X_BehavioralAnnotations', ...
        'sid_X_RawFrames', 'metadata', 'bdata', 'cdata', 'imagingData_raw', 'imagingData_reg'});

checkListTable = [table(expList.expID, 'VariableNames', {'ExpID'}), checkListTable];

writetable(checkListTable, fullfile(saveDir, 'preFT_expDirFiles.csv'))