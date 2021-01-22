

%% PLOT BACKGROUND ROI THROUGHOUT EACH EXPERIMENT
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\singleTrialAcq_pre_ficTrac';
roiFiles = dir(fullfile(parentDir, '*roiData.mat'));

for iExp = 1:numel(roiFiles)
    currFile = roiFiles(iExp).name;
    currExpID = currFile(1:10);
    disp(currExpID);

    load(fullfile(parentDir, currFile), 'roiData');
    backgroundRois = roiData(strcmp(roiData.roiName, 'Background'), :);
    figure(iExp); clf; 
    plot(smoothdata(cell2mat(backgroundRois.rawFl), 'gaussian', 3))
    title(currExpID)
    
end

%% Compare number of ROIdefs files to number of imgData block files

% Load exp list
expList = load_expList();

imgDataCounts = [];
roiDefCounts = [];
for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID);
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    parentDir = find_parent_dir(currExpID);
    
    imgDataFiles = dir(fullfile(parentDir, expDirName, 'imagingData_reg_block*.mat'));
    roiDefFiles = dir(fullfile(parentDir, expDirName, 'roiDefs_block*.mat'));
    
    if ~isempty(imgDataFiles)
        imgDataCounts(iExp) = numel(imgDataFiles);
    else
        imgDataCounts(iExp) = 0;
    end
    if ~isempty(roiDefFiles)
        roiDefCounts(iExp) = numel(roiDefFiles);
    else
        roiDefCounts(iExp) = 0;
    end
end

countTable = expList;
countTable.imgDataCounts = imgDataCounts';
countTable.roiDefCounts = roiDefCounts';
countDiff = imgDataCounts - roiDefCounts;
countTable.countDiff = countDiff';

%% Get a list of the ROI names for each experiment

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';
roiFiles = dir(fullfile(parentDir, '*roiData.mat'));

trialRoiNameTable = [];
for iExp = 1:numel(roiFiles)
    currFile = roiFiles(iExp).name;
    currExpID = currFile(1:10);
    disp(currExpID);

    load(fullfile(parentDir, currFile), 'roiData');
    trialRoiNameTable = [trialRoiNameTable; unique(roiData(:, 1:3))];
end

expRoiNameTable = unique(trialRoiNameTable(:, [1 3]));

    %% Check that each experiment has a 'background' ROI
    expList = unique(trialRoiNameTable(:, 1));
    disp(size(expList));
    disp(size(expRoiNameTable(strcmp(expRoiNameTable.roiName, 'Background'), :)))
    missingBGs = outerjoin(expList, expRoiNameTable(strcmp(expRoiNameTable.roiName, ...
            'Background'), :));

    %% Identify typos or other errors in ROI names
    uniqueNames = unique(expRoiNameTable.roiName);
    
    errorNames = expRoiNameTable(~cellfun(@isempty, ...
            regexp(expRoiNameTable.roiName, 'BackgroundDim')), :);

    %% Fix errors
    expID = '20200318-1';
    errorName = 'TypeD_2';
    correctName = 'TypeD-R-2';
    
    roiFileName = [expID, '_roiData.mat'];
    load(fullfile(parentDir, roiFileName), 'roiData');
    
    for iRow = 1:size(roiData, 1)
       if strcmp(errorName, roiData.roiName{iRow})
          roiData.roiName{iRow} = correctName; 
       end
    end
    
    save(fullfile(parentDir, roiFileName), 'roiData');
    clear roiFileName;
    
    



%%



nVolumes = 68;

roiTypes = unique(roiData(:, {'roiName'}));
% roiTypes = roiTypes(4, :);
for iType = 1:size(roiTypes, 1)
    currType = roiTypes.roiName{iType};
    currRoiData = innerjoin(roiData, roiTypes(iType, :));
    flMat = cell2mat(currRoiData.rawFl);
    flMat = reshape(flMat, nVolumes, []);
    figure(iType); clf
    imagesc(flMat');
    title(currType);
    figure(iType+10); clf
    plot(flMat);
    title(currType);
end

%%

cm = jet(size(test, 2));
figure(11);clf;hold on;
for i=1:size(test,2)
   plot(test(:, i), 'color', cm(i,:)) 
end

%%

roiTypes = unique(roiData(:, {'roiName'}));
for iType = 1:size(roiTypes, 1)
    currType = roiTypes.roiName{iType};
    currRoiData = innerjoin(roiData, roiTypes(iType, :));
    flMat = cell2mat(currRoiData.dffData');
    figure(iType); clf
    imagesc(flMat');
    title(currType);
end


%%

testData = currRoiData(1:20,:);
newTable = [];
for iRow = 1:size(testData, 1)
    currRow = testData(iRow, :);
    copyFields = currRow(:, [1:3, 5:7]);
    subROIs = currRow.subROIs{:};
    for iSub = 1:numel(subROIs)
       subRoiFields = struct2table(subROIs(iSub), 'asarray', 1);
       newRow = [copyFields, subRoiFields(:, 2:end)];
       newTable = [newTable; newRow];
    end
end

%%
