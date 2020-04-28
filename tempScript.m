
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\gaplessAcq';

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
allExpOdorEvents = odorEvent();
roiDataTable = [];
ftDataTable = [];
for iExp = 1:numel(gaplessAcqExpList)
    
    currExpID = gaplessAcqExpList{iExp};
    disp(currExpID);
    
    % trialMetadata files
    currTrialMd = readtable(fullfile(parentDir, [currExpID, '_trialMetadata.csv']));
    allTrialMd = [allTrialMd; currTrialMd];
    
    % Odor event data
    odorEventFile = fullfile(parentDir, [currExpID, '_event_data_odor.csv']);
    if exist(odorEventFile, 'file')
        allExpOdorEvents = allExpOdorEvents.load_csv(parentDir, 'fileNamePrefix', currExpID);
    end
    
    % ROI data
    load(fullfile(parentDir, [currExpID, '_roiData.mat']), 'roiData');
    for iTrial = 1:numel(roiData)
        for iROI = 1:numel(roiData{iTrial})
            newRow = table({currExpID}, 'VariableNames', {'expID'});
            newRow.trialNum = iTrial;
            newRow.roiName = {roiData{iTrial}(iROI).name};
            newRow.rawFl = {roiData{iTrial}(iROI).rawFl};
            newRow.expDffData = {roiData{iTrial}(iROI).expDffData};
            roiDataTable = [roiDataTable; newRow];
        end
    end% iTrial
    
    % FicTrac data
    load(fullfile(parentDir, [currExpID, '_ficTracData.mat']), 'ftData');    
    newData = [table(repmat({currExpID}, numel(ftData), 1), 'VariableNames', {'expID'}), ...
            struct2table(ftData, 'AsArray', 1)];
    ftDataTable = [ftDataTable; newData];
    
end%iExp


%%

combinedMetadata = innerjoin(gaplessAcqExpMd, allTrialMd);
oneNoteMetadata = readtable(fullfile(parentDir, 'OneNoteMetadata.txt'));

allMetadata = innerjoin(combinedMetadata, oneNoteMetadata);


test = innerjoin(roiDataTable, allMetadata(:, 1:end-2));



test2 = test( ...
        strcmp(test.roiName, 'ANT-R') & ... 
        test.starvationTimeHrs >= 24 & ...
        strcmp(test.bedtime, '4pm') ...
        , :);


%% SPLIT IMAGING DATA FILES INTO BLOCKS

parentDir = 'B:\ImagingData\';
expDirs = dir(fullfile(parentDir, '2*'));
expDirs = expDirs([expDirs.isdir]);
expDirs = expDirs(11:end); % Already did the first experiment

for iExp = 1:numel(expDirs)
    sidDir = dir(fullfile(expDirs(iExp).folder, expDirs(iExp).name, 's*'));
    sessionDataFile = dir(fullfile(sidDir.folder, sidDir.name, 'rigid*'));
    m = matfile(fullfile(sessionDataFile.folder, sessionDataFile.name));
    
    sz = size(m, 'wholeSession');
    
    if sz(5) <= 100
        tic
        disp('Loading entire session file...')
        load(fullfile(sessionDataFile.folder, sessionDataFile.name), 'wholeSession')
        disp(['Entire session loaded in ', num2str(toc, 3), ' sec'])
    end
    
    blockStartTrials = 1:20:sz(5);
    nBlocks = numel(blockStartTrials);
    for iBlock = 1:nBlocks
        disp([expDirs(iExp).name, ', block ', num2str(iBlock), ' of ', num2str(nBlocks)])
        if iBlock < nBlocks
            currBlockTrials = blockStartTrials(iBlock):(blockStartTrials(iBlock + 1) - 1);
        else
            currBlockTrials = blockStartTrials(iBlock):sz(5);
        end
        if sz(5) <= 100
           imgData = wholeSession(:, :, :, :, currBlockTrials); 
        else
            tic
            disp('Loading single block from matfile...')
            imgData = m.wholeSession(1:sz(1), 1:sz(2), 1:sz(3), 1:sz(4), currBlockTrials);
            disp(['Block loaded in ', num2str(toc, 3), ' sec'])
        end
        refImages = multi_mean(imgData, [4, 5]);
        
        tic
        save(fullfile(sessionDataFile.folder, ['imagingData_reg_block_', pad(num2str(iBlock), 2, ...
                'left', '0'), '.mat']), 'imgData');
        disp(['Block saved to disk in ', num2str(toc, 3), ' sec'])
        
        save(fullfile(sessionDataFile.folder, ['refImages_reg_block_', pad(num2str(iBlock), 2, ...
                'left', '0'), '.mat']), 'refImages');
            
        clear imgData
        
    end%iBlock
    
    clear wholeSession m
    
end%iExp


%% % CHECK NAMES OF FILES IN EACH EXPERIMENT DIRECTORY
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





