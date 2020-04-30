
%% PRE-FICTRAC EXPERIMENTS

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

%% SINGLE TRIAL ACQUISITION EXPERIMENTS

saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2018'; 
expList = readtable(fullfile(saveDir, 'oldExpList_singleTrialAcq.csv'), 'delimiter', ',');

allFilesList = [];
dirInfoTable = []; dirCheckList = [];
for iExp = 1:size(expList, 1)
   
    currExpID = expList.expID{iExp};
    disp(currExpID);
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    
    % Find experiment directory
    if isempty(dir(fullfile(parentDir, expDirName)))
        oldExpParentDirs = dir(fullfile(parentDir, '2018 *'));
        for iDir = 1:numel(oldExpParentDirs)
            expParentDir = fullfile(oldExpParentDirs(iDir).folder, oldExpParentDirs(iDir).name);
            if ~isempty(dir(fullfile(expParentDir, expDirName)))
                break
            end
            expParentDir = 'F:\ImagingData';
        end
    else
        expParentDir = fullfile(parentDir);
    end
    sidDir = dir(fullfile(expParentDir, expDirName, 'sid_*'));
    sidDir = sidDir([sidDir.isdir]);
    if numel(sidDir) == 1
        expDir = fullfile(expParentDir, expDirName, sidDir.name); % If there's only one sid directory
    elseif numel(sidDir) > 1
        expDir = fullfile(expParentDir, expDirName, 'sid_master'); % If I've already manually combined sids
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
    
    dirCheckList(iExp, 1) = exist(fullfile(expDir, 'allTrials.csv'), 'file');
    dirCheckList(iExp, 2) = exist(fullfile(expDir, 'Annotations.mat'), 'file');
    dirCheckList(iExp, 3) = numel(dir(fullfile(expDir, 'sid_0_Annotations.mat'))) > 0;
    dirCheckList(iExp, 4) = exist(fullfile(expDir, 'analysisMetadata.mat'), 'file');
    dirCheckList(iExp, 5) = numel(dir(fullfile(expDir, 'ROI_metadata*.mat'))) > 0;
    dirCheckList(iExp, 6) = numel(dir(fullfile(expDir, 'ROImetadata.mat'))) > 0;
    dirCheckList(iExp, 7) = numel(dir(fullfile(expDir, 'refImages_Reg.mat'))) > 0;
    dirCheckList(iExp, 8) = numel(dir(fullfile(expDir, 'refImages_s*.mat'))) > 0;
    dirCheckList(iExp, 9) = numel(dir(fullfile(expDir, 'sid*refImages.mat'))) > 0;
    dirCheckList(iExp, 10) = numel(dir(fullfile(expDir, '*volAvgSessionData.mat'))) > 0;
    dirCheckList(iExp, 11) = numel(dir(fullfile(expDir, '*frameCountLog.mat'))) > 0;
    dirCheckList(iExp, 12) = numel(dir(fullfile(expDir, 'sid_*_optic_flow_data.mat'))) > 0;
    dirCheckList(iExp, 13) = numel(dir(fullfile(expDir, '*flow_data_norm.mat'))) > 0;
    dirCheckList(iExp, 14) = numel(dir(fullfile(expDir, 'skipTrials.mat'))) > 0;
    dirCheckList(iExp, 15) = numel(dir(fullfile(expDir, 'ROI_Data_Avg.mat'))) > 0;
    dirCheckList(iExp, 16) = numel(dir(fullfile(expDir, 'annotationTypes.mat'))) > 0;
    dirCheckList(iExp, 17) = numel(dir(fullfile(expDir, 'annotationParams.mat'))) > 0;
    dirCheckList(iExp, 18) = numel(dir(fullfile(expDir, 'imagingData_reg_block*.mat'))) > 0;
    dirCheckList(iExp, 19) = numel(dir(fullfile(expDir, 'imgMetadata.mat'))) > 0;
    dirCheckList(iExp, 20) = numel(dir(fullfile(expDir, 'refImages_reg_block*.mat'))) > 0;
    dirCheckList(iExp, 21) = numel(dir(fullfile(expDir, 'refImages_sid_*_sessionFile.mat'))) > 0;
    dirCheckList(iExp, 22) = numel(dir(fullfile(expDir, 'rigid_wholeSession.mat'))) > 0;
    dirCheckList(iExp, 23) = numel(dir(fullfile(expDir, 'sid*allTrials_frameCountLog.mat'))) > 0;
    dirCheckList(iExp, 24) = numel(dir(fullfile(expDir, 'sid*sessionFile.mat'))) > 0;  
    dirCheckList(iExp, 25) = numel(dir(fullfile(expDir, 'metadata*.mat'))) > 0;  
end

dirCheckList = logical(dirCheckList);
checkListTable = array2table(dirCheckList, 'VariableNames', {'allTrials', 'Annotations', ...
        'sid_0_Annotations', 'analysisMetadata', 'ROI_metadata', 'ROImetadata', 'refImages_Reg', ...
        'refImages_sid', 'sid_X_refImages', 'volAvgSessionData', 'frameCountLog', ...
        'sid_X_optic_flow_data', 'flow_data_norm', 'skipTrials', 'ROI_Data_Avg', 'annotationTypes', ...
        'annotationParams', 'imagingData_reg_block_X', 'imgMetadata', 'refImages_reg_block_X', ...
        'refImages_sid_X_sessionFile', 'rigid_wholeSession', 'allTrials_frameCountLog', ...
        'sid_X_sessionFile', 'metadata'});

checkListTable = [table(expList.expID, 'VariableNames', {'ExpID'}), checkListTable];

% writetable(checkListTable, fullfile(saveDir, 'expDirFiles.csv'))

%% CHECK SID COUNTS

parentDir = 'F:\ImagingData';
expDirs = dir(fullfile(parentDir, '20*'));
expDirs = expDirs([expDirs.isdir]);

for iExp = 1:numel(expDirs)
    expDirContents = dir(fullfile(parentDir, expDirs(iExp).name));
    expDirContents = {expDirContents.name}';
    expDirContents = expDirContents(~cellfun(@(x) (strcmp(x, '.') | strcmp(x, '..')), ...
            expDirContents));
    disp(' ')
    disp(expDirs(iExp).name)
    disp(expDirContents)
end







