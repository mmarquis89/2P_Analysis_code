
saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';

expID = '20200318-2';

expDir = dir(fullfile(parentDir, [expID, '*']));
sourceDataDir = fullfile(expDir.folder, expDir.name, 'ProcessedData');


% Load necessary files
load(fullfile(sourceDataDir, 'analysis_data.mat'), 'mD');
load(fullfile(sourceDataDir, 'flowMags.mat'), 'meanFlowMags');



% ----------- Create variables for new files --------------

% Experiment metadata
expMetadata = table({expID}, {mD(1).expMetadata.expName}, mD(1).daqSampRate, mD(1).panelsDisplayRate, ...
        mD(1).volumeRate, mD(1).nDaqSamples, mD(1).nPanelsFrames, mD(1).nVolumes, ...
        mD(1).trialDuration, mD(1).nPlanes, numel(mD));
expMetadata.Properties.VariableNames = {'expID', 'expName', 'daqSampRate', ...
        'panelsDisplayRate', 'volumeRate', 'nDaqSamples', 'nPanelsFrames', 'nVolumes', ...
        'trialDuration', 'nPlanes', 'nTrials'};

% Trial metadata
mdTable = struct2table(mD);
trialMetadata = mdTable(:, {'trialNum', 'usingPanels', 'using2P'}); 

% Panels metadata
panelsMetadata = struct();
for iTrial = 1:expMetadata.nTrials
    panelsMetadata(iTrial).trialNum = trialMetadata.trialNum(iTrial);
    panelsMetadata(iTrial).panelsMode = mD(iTrial).expMetadata.panelsMode;
    panelsMetadata(iTrial).pattern = mD(iTrial).expMetadata.pattern; 
    panelsMetadata(iTrial).xPosFunc = mD(iTrial).expMetadata.xDimPosFun.func;
    panelsMetadata(iTrial).yPosFunc = mD(iTrial).expMetadata.yDimPosFun.func;
end
   
% FicTrac data 
ftData = struct;
for iTrial = 1:expMetadata.nTrials
    ftData(iTrial).trialNum = trialMetadata.trialNum(iTrial);
    ftData(iTrial).intX = mD(iTrial).ftData.intX;
    ftData(iTrial).intY = mD(iTrial).ftData.intY;
    ftData(iTrial).intHD = mD(iTrial).ftData.intHD;
    ftData(iTrial).moveSpeed = mD(iTrial).ftData.moveSpeed;
    ftData(iTrial).intFwMove = mD(iTrial).ftData.intFwMove;
    ftData(iTrial).intSideMove = mD(iTrial).ftData.intSideMove;
    ftData(iTrial).yawSpeed = mD(iTrial).ftData.yawSpeed;
    ftData(iTrial).fwSpeed = mD(iTrial).ftData.fwSpeed;
    ftData(iTrial).sideSpeed = mD(iTrial).ftData.sideSpeed;
    ftData(iTrial).badFtTrials = nan;
end 
    
% ROI defs
roiDefs = {};
for iTrial = 1:expMetadata.nTrials
    try
        roiDefs{iTrial} = rmfield(mD(iTrial).roiData, {'color', 'rawData', 'dffData', 'zscoreData', ...
                'expDffData'});
    catch
        roiDefs{iTrial} = rmfield(mD(iTrial).roiData, {'color', 'rawData', 'dffData', 'zscoreData'});
    end
end
    
% ----------- Save data in group analysis directory --------------
writetable(expMetadata, fullfile(saveDir, ['ExpMetadata_', expID, '.csv'])); 
writetable(trialMetadata, fullfile(saveDir, ['TrialMetadata_', expID, '.csv'])); 
save(fullfile(saveDir, ['PanelsMetadata_', expID, '.mat']), 'panelsMetadata'); 
save(fullfile(saveDir, ['FicTracData_', expID, '.mat']), 'ftData'); 
save(fullfile(saveDir, ['ROI_defs_', expID, '.mat']), 'roiDefs'); 










