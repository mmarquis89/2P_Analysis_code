% Need the following fields:
% 
%   volumeRate
%   nPlanes
%   nTrials
%   trialDuration
%   nVolumes
%   nFrames
% refImg
%
%

parentDir = 'F:\ImagingData'; 
expList = readtable(fullfile('D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData', ...
        'oldExpList_singleTrialAcq.csv'), 'delimiter', ',');
expList = expList(1:10, :);

for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID);
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    
    dataDir = dir(fullfile(parentDir, expDirName, 'sid_*'));
    if numel(dataDir) > 1
        error('Multiple sid folders detected')
    else
        dataDir = fullfile(dataDir(1).folder, dataDir(1).name);
    end
    
    aD = [];
    aD.expDate = expDirName;
    
%     % DAQ metadata files
%     mdFiles = dir(fullfile(dataDir, 'metadata*.mat'));
%     aD.nTrials = numel(mdFiles);
%     load(fullfile(mdFiles(1).folder, mdFiles(1).name), 'metaData');
%     aD.trialDuration = metaData.trialDuration;
    
    % Imaging data (just to check size)
    m = matfile(fullfile(dataDir, 'rigid_wholeSession.mat'));
    try
        sz = size(m, 'regProduct');
    catch
       sz = size(m, 'wholeSession'); 
    end
    aD.nPlanes = sz(3);
    aD.nVolumes = sz(4);
    
    % SI metadata
    try
        load(fullfile(dataDir, 'imgMetadata.mat'));
        siInfo = imgMetadata.scanimageInfo;
    catch
        siInfo = m.scanimageInfo;
    end
    volRateInfo = regexp(siInfo, '(?<=scanimage.SI.hRoiManager.scanVolumeRate = ).*', 'match');
    volRateInfo = volRateInfo{~cellfun(@isempty, volRateInfo)};
    aD.volumeRate = str2double(volRateInfo{:});
    clear m;
    aD.trialDuration = round(aD.nVolumes / aD.volumeRate);
    
    % Annotation data (just for frame count)
    load(fullfile(dataDir, 'Annotations.mat'), 'trialAnnotations')
    aD.nFrames = mode(cellfun(@(x) size(x, 1), trialAnnotations));
    
    % Reference images
    refImgFile = dir(fullfile(dataDir, 'refImages*.mat'));
    if numel(refImgFile) > 1
        error('Multiple refImgFiles detected')
    else
        refImgFile = fullfile(refImgFile(1).folder, refImgFile(1).name);
    end
    load(refImgFile, 'refImages')
    aD.refImg = refImages;

    
    analysisMetadata = aD;
    save(fullfile(dataDir, 'analysisMetadata.mat'), 'analysisMetadata');
end





