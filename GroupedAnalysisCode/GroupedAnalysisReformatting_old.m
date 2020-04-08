
saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2019 Mar-Apr\2019_04_01_exp_2\sid_0';
expName = 'D-ANT_6s';


% Load analysisMetadata file
load(fullfile(parentDir, 'analysisMetadata.mat', analysisMetadata));
aD = analysisMetadata;

FRAME_RATE = 25; % Behavior video (and therefore FicTrac) frame rate was a constant

% ----- CREATE EXP METADATA TABLE -----

expID = regexprep(aD.expDate, {'_', 'exp'}, {'', '-'});
expMd = table({expID}, 'VariableNames', {'expID'});
expMd.expName = expName;
expMd.daqSampRate = 4000;        % (this is after 10x downsampling at the time of the experiment)
expMd.panelsDisplayRate = 50;    % Didn't have panels but just keeping the number consistent
expMd.volumeRate = aD.volumeRate;
expMd.nPlanes = aD.nPlanes; 
expMd.nTrials = numel(aD.blockData); % Redefining "trial" as an individual continuous acquisition

% ----- CREATE TRIAL METADATA TABLE -----
trialMd = [];
for iTrial = 1:expMd.nTrials
    newRow = table(iTrial, 'VariableNames', {'trialNum'});
    newRow.trialDuration = aD.trialDuration * aD.blockData(iTrial).nTrials;
    newRow.nVolumes = aD.nVolumes * aD.blockData(iTrial).nTrials;
    newRow.nDaqSamples = size(aD.blockData.outputData, 1);
    newRow.nPanelsFrames = expMd.panelsDisplayRate * newRow.trialDuration;
    newRow.usingOptoStim = 0;
    newRow.usingPanels = 0;
    newRow.using2P = 1;
    trialMd = [trialMd; newRow];
end

% ----- CREATE FICTRAC DATA STRUCTURE -----
nBlockTrials = [aD.blockData.nTrials];
blockStartTrials = 1 + [0, cumsum(nBlockTrials(1:end-1))];
blockEndTrials = cumsum(nBlockTrials);   
ftData = struct();
for iBlock = 1:numel(nBlockTrials)
    currBlockTrials = blockStartTrials(iBlock):blockEndTrials(iBlock);
    ftData(iBlock).trialNum = iBlock;
    ftData(iBlock).intX = as_vector(squeeze(aD.ftData.intXY(:, 1, currBlockTrials)));
    ftData(iBlock).intY = as_vector(squeeze(aD.ftData.intXY(:, 2, currBlockTrials)));
    ftData(iBlock).intHD = as_vector(aD.ftData.intHD(:, currBlockTrials));
    ftData(iBlock).moveSpeed = as_vector(aD.ftData.moveSpeed(:, currBlockTrials));
    ftData(iBlock).intFwMove = as_vector(aD.ftData.intForwardMove(:, currBlockTrials));
    ftData(iBlock).intSideMove = as_vector(aD.ftData.intSideMove(:, currBlockTrials));
    ftData(iBlock).yawSpeed = as_vector(aD.ftData.yawSpeed(:, currBlockTrials));
    ftData(iBlock).fwSpeed = as_vector(aD.ftData.fwSpeed(:, currBlockTrials));
    ftData(iBlock).sideSpeed = as_vector(aD.ftData.sideSpeed(:, currBlockTrials));
    currBlockFrameCount = numel(ftData(iBlock).intX);
    ftData(iBlock).frameTimes = (1:currBlockFrameCount) / FRAME_RATE;
    ftData(iBlock).badVidFrames = isnan(ftData(iBlock).intX);
    ftData(iBlock).meanFlow = as_vector(aD.flowArr(:, currBlockTrials));
end


% ----- 








