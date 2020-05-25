
%% Blocks of processed imaging data

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017\2017_11_27_exp_2\sid_2';

imgDataFiles = dir(fullfile(parentDir, 'imagingData_reg_block*.mat'));
medVals = [];
minVals = [];
meanVals = [];
for iFile = 1:numel(imgDataFiles)
    disp(iFile)
    load(fullfile(parentDir, imgDataFiles(iFile).name), 'imgData'); % --> [y, x, plane, vol, trial]
    sz = size(imgData);
    imgDataPerm = permute(reshape(imgData, sz(1), sz(2), sz(3), []), [4 1 2 3]); % --> [expVol, y, x, plane]
    medVals = [medVals, median(reshape(imgDataPerm, size(imgDataPerm, 1), []), 2)'];
    minVals = [minVals, min(reshape(imgDataPerm, size(imgDataPerm, 1), []), [], 2)'];
    meanVals = [meanVals, mean(reshape(imgDataPerm, size(imgDataPerm, 1), []), 2)'];
end


%% Cdata .mat files

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017\2017_11_27_exp_2\sid_2';

planeNum = 6;

imgDataFiles = dir(fullfile(parentDir, 'cdata*.mat'));
medVals = [];
minVals = [];
meanVals = [];
for iFile = 1:numel(imgDataFiles)
    disp(iFile)
    load(fullfile(parentDir, imgDataFiles(iFile).name), 'tifData'); % --> [y, x, plane, vol]
    tifDataPerm = permute(tifData, [4 1 2 3]);                      % --> [vol, y, x, plane]
    medVals = [medVals, median(reshape(tifDataPerm, size(tifDataPerm, 1), []), 2)'];
    minVals = [minVals, min(reshape(tifDataPerm, size(tifDataPerm, 1), []), [], 2)'];
    meanVals = [meanVals, mean(reshape(tifDataPerm, size(tifDataPerm, 1), []), 2)'];
end

%%

sz = size(wholeSession); 
wholeSession = reshape(wholeSession, sz(1), sz(2), sz(3), []);  % --> [y, x, plane, expVol]
wholeSession = permute(wholeSession, [4 1 2 3]);                % --> [expVol, y, x, plane]
wholeSession = reshape(wholeSession, size(wholeSession, 1), []);% --> [expVol, px]
medSessionVals = median(wholeSession, 2);
minSessionVals = min(wholeSession, [], 2);
meanSessionVals = mean(wholeSession, 2);


%%

sz = size(wholeSession);
wholeSession = reshape(wholeSession, sz(1), sz(2), sz(3), []);  % --> [y, x, plane, expVol]
wholeSession = permute(wholeSession, [4 1 2 3]);                % --> [expVol, y, x, plane]
wholeSession = reshape(wholeSession, size(wholeSession, 1), []);% --> [expVol, px]
medRegVals = median(wholeSession, 2);
minRegVals = min(wholeSession, [], 2);
meanRegVals = mean(wholeSession, 2);

%%

xx = 1:129:size(meanVals, 2);

figure(1); clf; hold on
plot(minVals)
plot(medVals)
plot(smoothdata(meanVals, 'gaussian', 7))
plot(xx, ones(size(xx)) * 20, 'o', 'color', 'm')

% figure(2); clf; hold on
% plot(minSessionVals);
% plot(medSessionVals);
% plot(smoothdata(meanSessionVals, 'gaussian', 11));
% plot(xx, ones(size(xx)) * 50, 'o')
% 
% figure(3); clf; hold on
% plot(minRegVals);
% plot(medRegVals);
% plot(smoothdata(meanRegVals, 'gaussian', 11));
% plot(xx, ones(size(xx)) * 50, 'o')

%%

sz = size(imgData);
rsData = reshape(imgData, sz(1), sz(2), sz(3), []); % [y, x, plane, vol]
permData = permute(rsData, [4 1 2 3]); % [vol, y, x, plane]
sz = size(permData);
medPxByVol = median(reshape(permData, sz(1), []), 2); % [vol]
sz = size(rsData);
medPxRep = permute(repmat(medPxByVol, 1, sz(1), sz(2), sz(3)), [2 3 4 1]); % [y, x, plane, volume]

baseSubData = rsData - medPxRep;
meanVals = squeeze(multi_mean(baseSubData, [1 2]));        
sz = size(baseSubData);
medVals = squeeze(median(reshape(permute(baseSubData, [3 4 1 2]), sz(3), sz(4), []), 3)); % --> [plane, vol] 
minVals = squeeze(min(reshape(permute(baseSubData, [3 4 1 2]), sz(3), sz(4), []), [], 3));% --> [plane, vol] 


%% PLOT ROI DATA TO CHECK FOR OFFSETS AT TRIAL BOUNDS

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\testing';

currFileName = 'newFormat.mat'; %'gaplessAcq.mat'; %'preGaplessAcq.mat'; %  

smWin = 7;

load(fullfile(parentDir, currFileName), 'allExpROIData');
expIDs = unique(allExpROIData.expID);
for iExp = 1:numel(expIDs)
   currData = allExpROIData(strcmp(allExpROIData.expID, expIDs{iExp}), :); 
   
   % Create figure
   f = figure(iExp); clf; hold on
   f.Color = [1 1 1];
   f.Position = [50, 250, 1500, 700];
   
   for iROI = 1:size(currData, 1)
       flData = currData.rawFl{iROI}(:);
       flData(currData.trialBounds{1}) = nan;
       plot(smoothdata(flData, 'gaussian', smWin, 'omitnan'));
   end
   legend(currData.roiName, 'location', 'best', 'autoupdate', 'off')
   yL = ylim();
   for iTrial = 1:numel(currData.trialBounds{1})
       startVol = currData.trialBounds{1}(iTrial);
       plot([startVol, startVol], yL, '--', 'color', 'm')
   end
   
   
end

%% CHECK MEDIAN VALUES FOR PRE-ROI EXPERIMENTS

parentDir = 'F:\ImagingData';
expList = load_expList('groupName', 'singleTrialAcq_preROIs');

for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID)
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    
    blockFiles = dir(fullfile(parentDir, expDirName, 'sid_0', 'imagingData_reg_block*.mat'));
    medVals = [];
    for iFile = 1:numel(blockFiles)
        disp(iFile)
        load(fullfile(blockFiles(iFile).folder, blockFiles(iFile).name), 'imgData');
        sz = size(imgData);
        nVolumes = sz(4);
        rsData = reshape(imgData, sz(1), sz(2), sz(3), []);
        sz = size(rsData);
        medVals = [medVals, median(reshape(permute(rsData, [4 1 2 3]), sz(4), []), 2)'];        
    end
    expList.medVals{iExp} = medVals;
    expList.nVolumes(iExp) = nVolumes;
end

%% CHECK MEDIAN VALUES FOR PRE-FICTRAC 2018 EXPERIMENTS

expList = load_expList('groupName', 'singleTrialAcq_preFicTrac');
expList = expList(17:end, :)
    
    
for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID)
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    sidDir = dir(fullfile(parentDir, expDirName, 'sid_*'));
    if numel(sidDir) > 1
        blockFiles = dir(fullfile(parentDir, expDirName, 'sid_master', 'imagingData_reg_block*.mat'));
    else
        blockFiles = dir(fullfile(parentDir, expDirName, sidDir(1).name, 'imagingData_reg_block*.mat'));
    end
    
    medVals = [];
    for iFile = 1:numel(blockFiles)
        disp(iFile)
        try
        load(fullfile(blockFiles(iFile).folder, blockFiles(iFile).name), 'imgData');
        catch
            
        end
        sz = size(imgData);
        nVolumes = sz(4);
        rsData = reshape(imgData, sz(1), sz(2), sz(3), []);
        sz = size(rsData);
        medVals = [medVals, median(reshape(permute(rsData, [4 1 2 3]), sz(4), []), 2)'];        
    end
    expList.medVals{iExp} = medVals;
    expList.nVolumes(iExp) = nVolumes;
end
%% CHECK MEDIAN VALUES FOR PRE-FICTRAC 2017 EXPERIMENTS

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017';
expList = load_expList('groupName', 'singleTrialAcq_preFicTrac');
expList = expList(3:16, :)
    
    
for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID)
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    sidDir = dir(fullfile(parentDir, expDirName, 'sid_*'));
    if numel(sidDir) > 1
        blockFiles = dir(fullfile(parentDir, expDirName, 'sid_master', 'imagingData_reg_block*.mat'));
    else
        blockFiles = dir(fullfile(parentDir, expDirName, sidDir(1).name, 'imagingData_reg_block*.mat'));
    end
    
    medVals = [];
    for iFile = 1:numel(blockFiles)
        disp(iFile)
        try
        load(fullfile(blockFiles(iFile).folder, blockFiles(iFile).name), 'imgData');
        catch
            
        end
        sz = size(imgData);
        nVolumes = sz(4);
        rsData = reshape(imgData, sz(1), sz(2), sz(3), []);
        sz = size(rsData);
        medVals = [medVals, median(reshape(permute(rsData, [4 1 2 3]), sz(4), []), 2)'];        
    end
    expList.medVals{iExp} = medVals;
    expList.nVolumes(iExp) = nVolumes;
end
%%
for iExp = 1:size(expList, 1)
    
   % Create figure
   f = figure(iExp); clf; hold on
   f.Color = [1 1 1];
   f.Position = [50, 250, 1500, 700];
   plot(smoothdata(expList.medVals{iExp}, 'gaussian', 5));
   yL = ylim();
   trialBounds = 1:expList.nVolumes(iExp):numel(expList.medVals{iExp});
   for iTrial = 1:numel(trialBounds)
       startVol = trialBounds(iTrial);
       plot([startVol, startVol], yL, '--', 'color', 'm')
   end
end

%% Load ROI data to compare offset across trials


% PRE-GAPLESS ACQUISITION
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\singleTrialAcq_preROIs'; 
roiFiles = dir(fullfile(parentDir, '*_roiData.mat'));

allExpROIData = [];
for iExp = 1:numel(roiFiles)
    currFileName = roiFiles(iExp).name;
    currExpID = currFileName(1:10);
    disp(currExpID)
    newRow = table({currExpID}, 'VariableNames', {'expID'});
%     % Load trial metadata
%     trialMd = readtable(fullfile(saveDir, [currExpID, '_trialMetadata.csv']), 'delimiter', ',');
    
    % Load ROI data
    load(fullfile(parentDir, currFileName), 'roiData');
    roiNames = unique(roiData.roiName);
    nROIs = numel(roiNames);
    for iROI = 1:nROIs
        newRow.roiName = roiNames(iROI);
        currROIData = roiData(strcmp(roiNames{iROI}, roiData.roiName), :);
        fullExpFlData = [];
        trialBounds = 1;
        for iTrial = 1:size(currROIData, 1)
            fullExpFlData = [fullExpFlData, currROIData.rawFl{iTrial}];
            trialBounds(end + 1) = trialBounds(end) + numel(currROIData.rawFl{iTrial});
        end
        newRow.rawFl = {fullExpFlData};
        newRow.trialBounds = {trialBounds(1:end-1)};
        allExpROIData = [allExpROIData; newRow];
    end
end


% OLD FORMAT, GAPLESS ACQUISITION
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\gaplessAcq'; 
roiFiles = dir(fullfile('D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData', '*_roiData.mat'));

allExpROIData = [];
for iExp = 1:numel(roiFiles)
    currFileName = roiFiles(iExp).name;
    currExpID = currFileName(1:10);
    disp(currExpID)
    newRow = table({currExpID}, 'VariableNames', {'expID'});
    
    % Load trial metadata
    trialMd = readtable(fullfile(parentDir, [currExpID, '_trialMetadata.csv']), 'delimiter', ',');
    
    % Load ROI data
%     load(fullfile(parentDir, currFileName), 'roiData');
    load(fullfile('D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData', currFileName), 'roiData');
    roiData = innerjoin(roiData, trialMd(:, {'expID', 'trialNum', 'originalTrialCount'}));
    roiNames = unique(roiData.roiName);
    nROIs = numel(roiNames);
    for iROI = 1:nROIs
        newRow.roiName = roiNames(iROI);
        currROIData = roiData(strcmp(roiNames{iROI}, roiData.roiName), :);
        fullExpFlData = [];
        trialBounds = 1;
        for iTrial = 1:size(currROIData, 1)
            volCount = numel(currROIData.rawFl{iTrial}) / currROIData.originalTrialCount(iTrial);
            fullExpFlData = [fullExpFlData, currROIData.rawFl{iTrial}'];
            trialBounds = [trialBounds, trialBounds(end) + ...
                    (volCount:volCount:numel(currROIData.rawFl{iTrial}))];
        end
        newRow.rawFl = {fullExpFlData};
        newRow.trialBounds = {trialBounds(1:end-1)};
        allExpROIData = [allExpROIData; newRow];
    end
end

% NEW FORMAT
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\newData'; 
roiFiles = dir(fullfile(parentDir, '*_roiData.mat'));

allExpROIData = [];
for iExp = 1:numel(roiFiles)
    currFileName = roiFiles(iExp).name;
    currExpID = currFileName(1:10);
    disp(currExpID)
    newRow = table({currExpID}, 'VariableNames', {'expID'});
    
    % Load ROI data
    load(fullfile(parentDir, currFileName), 'roiData');
    roiNames = unique(roiData.roiName);
    nROIs = numel(roiNames);
    for iROI = 1:nROIs
        newRow.roiName = roiNames(iROI);
        currROIData = roiData(strcmp(roiNames{iROI}, roiData.roiName), :);
        fullExpFlData = [];
        trialBounds = 1;
        for iTrial = 1:size(currROIData, 1)
            fullExpFlData = [fullExpFlData, currROIData.rawFl{iTrial}];
            trialBounds(end + 1) = trialBounds(end) + numel(currROIData.rawFl{iTrial});
        end
        newRow.rawFl = {fullExpFlData};
        newRow.trialBounds = {trialBounds(1:end-1)};
        allExpROIData = [allExpROIData; newRow];
    end
end


