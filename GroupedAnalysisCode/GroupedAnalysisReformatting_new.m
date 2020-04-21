
saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';

expNum = 27;

try
    
expName = expNames.expDirs{expNum};
sourceDataDir = fullfile(parentDir, expName, 'ProcessedData');

% expDir = uigetdir(parentDir, 'Select an experiment directory');
% sourceDataDir = fullfile(expDir, 'ProcessedData');
% expID = regexp(sourceDataDir, '(?<=Imaging Data\\)........-.(?=_)', 'match', 'once');
% expID = '20200127-2';

% expDir = dir(fullfile(parentDir, [expID, '*']));
% sourceDataDir = fullfile(expDir.folder, expDir.name, 'ProcessedData');


% Load necessary files
load(fullfile(sourceDataDir, 'analysis_data.mat'), 'mD');
load(fullfile(sourceDataDir, 'flowMags.mat'), 'meanFlowMags');

expID = mD(1).expID;

% ----------- Create variables for new files --------------

% Experiment metadata
expMetadata = table({expID}, {mD(1).expMetadata.expName}, mD(1).daqSampRate, mD(1).panelsDisplayRate, ...
        mD(1).volumeRate, mD(1).nPlanes, numel(mD));
expMetadata.Properties.VariableNames = { ...
        'expID', ... 
        'expName', ...
        'daqSampRate', ...
        'panelsDisplayRate', ...
        'volumeRate', ...
        'nPlanes', ...
        'nTrials'};

% Trial metadata
mdTable = struct2table(mD);
trialMetadata = mdTable(:, {'trialNum', 'trialDuration', 'nVolumes', 'nDaqSamples', 'nPanelsFrames', ...
        'usingOptoStim', 'usingPanels', 'using2P'}); 
trialMetadata = [table(repmat({expID}, size(trialMetadata, 1), 1), 'VariableNames', {'expID'}), ...
        trialMetadata];

% Panels metadata
panelsMetadata = struct();
for iTrial = 1:expMetadata.nTrials
    panelsMetadata(iTrial).trialNum = trialMetadata.trialNum(iTrial);
    panelsMetadata(iTrial).panelsMode = mD(iTrial).expMetadata.panelsMode;
    panelsMetadata(iTrial).pattern = mD(iTrial).expMetadata.pattern; 
    panelsMetadata(iTrial).xPosFunc = mD(iTrial).expMetadata.xDimPosFun.func;
    panelsMetadata(iTrial).yPosFunc = mD(iTrial).expMetadata.yDimPosFun.func;
    panelsMetadata(iTrial).panelsPosX = mD(iTrial).panelsPosX;
    panelsMetadata(iTrial).panelsPosY = mD(iTrial).panelsPosY;
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
    ftData(iTrial).frameTimes = mD(iTrial).ftFrameTimes;
    ftData(iTrial).badFtTrials = nan;
    ftData(iTrial).meanFlow = meanFlowMags{iTrial};
end 
    
% ROI defs and extracted data
roiData = {};
for iTrial = 1:expMetadata.nTrials
    if ~isempty(mD(iTrial).roiData)
        roiData{iTrial} = rmfield(mD(iTrial).roiData, {'color', 'zscoreData'});
    else
        roiData{iTrial} = [];
    end
end

% Reference images 
refImages = zeros([size(mD(1).refImages), expMetadata.nTrials]);
for iTrial = 1:expMetadata.nTrials
    refImages(:, :, :, iTrial) = mD(iTrial).refImages; % --> [y, x, plane, trial]
end


% ----------- Save data in group analysis directory --------------
writetable(expMetadata, fullfile(saveDir, [expID, '_expMetadata.csv'])); 
writetable(trialMetadata, fullfile(saveDir, [expID, '_trialMetadata.csv'])); 
save(fullfile(saveDir, [expID, '_panelsMetadata.mat']), 'panelsMetadata'); 
save(fullfile(saveDir, [expID, '_ficTracData.mat']), 'ftData'); 
save(fullfile(saveDir, [expID, '_ROI_data.mat']), 'roiData'); 
save(fullfile(savedir, [expID, '_refImages.mat']), 'refImages');


% CREATE EVENT OBJECTS AND EXPORT DATA

% ----- Panels flash stim -----

if any(trialMetadata.usingPanels)
    panelsFlashStims = panelsFlashEvent();
    panelsTrials = find(trialMetadata.usingPanels);
    for iTrial = 1:numel(panelsTrials)
        currTrial = panelsTrials(iTrial);
        panelsPosX = panelsMetadata(currTrial).panelsPosX;
        
        % Use pattern file name to figure out what stim type(s) were shown, and the brightness level
        patternName = panelsMetadata(currTrial).pattern.userData.patternName;
        if contains(patternName, 'All-Panels-On')
            patternTypes = {'fullField'};
            patternPositions = max(panelsPosX);
            brightness = round(max(panelsPosX) * (100/14));
        elseif contains(patternName,  'Ydim-brightness') || contains(patternName, 'bright_bar')
            if numel(unique(panelsMetadata(currTrial).xPosFunc)) > 50
                % Skip trial if it's one of the very few swinging bar stims
                continue
            end 
            patternTypes = {'leftBar', 'centerBar', 'rightBar'};
            patternPositions = [19 43 67];
            brightness = round(max(panelsMetadata(currTrial).panelsPosY) * (100 /14));
        else
            error('Invalid pattern')
        end
        
        % Calculate event times for each stim type and add to event object
        for iType = 1:numel(patternTypes)
            stimFrames = panelsPosX == patternPositions(iType);
            if sum(stimFrames)
                stimOnsetFrames = find(diff(stimFrames) == 1) + 1;
                stimOffsetFrames = find(diff(stimFrames) == -1) + 1;
                frameTimes = (1:trialMetadata(iTrial).nPanelsFrames) / ...
                        expMetadata.panelsDisplayRate;
                stimOnsetTimes = frameTimes(stimOnsetFrames);
                stimOffsetTimes = frameTimes(stimOffsetFrames);
                
                trialNum = panelsMetadata(currTrial).trialNum;
                eventTimes = {[stimOnsetTimes', stimOffsetTimes']};
                
                panelsFlashStims = panelsFlashStims.append_data(expMetadata.expID{1}, trialNum, ...
                        eventTimes, brightness, patternTypes{iType});
            end
        end%iType
        
    end%iTrial
    
    % Export to .csv file
    panelsFlashStims.export_csv(saveDir, 'fileNamePrefix', expMetadata.expID{1})
    
    
end%if 

catch ME
    rethrow(ME);
end

%% ----- Odor stim -----

odorNames = {'ACV'};
concentrations = {'neat'};
flowRates = {12};
trialNums = {[]};

try
if any(trialMetadata.usingOptoStim)
    odorStims = odorEvent();
    odorTrials = find(trialMetadata.usingOptoStim);
    for iCond = 1:numel(trialNums)
        
        % If there are multiple types of stim delivery append them one at a time
        odorName = odorNames{iCond};
        concentration = concentrations{iCond};
        flowRate = flowRates{iCond};
        currTrialNums = trialNums{iCond};
        if isempty(currTrialNums)
            currTrialNums = trialMetadata.trialNum;
        end
        for iTrial = 1:numel(currTrialNums)
            if trialMetadata.usingOptoStim(trialMetadata.trialNum == currTrialNums(iTrial))
                
                stimTiming = mD(trialMetadata.trialNum == currTrialNums(iTrial)).optoStimTiming;
                trialDuration = trialMetadata.trialDuration(trialMetadata.trialNum == ...
                        currTrialNums(iTrial));
                odorStims = odorStims.append_shorthand(expMetadata.expID{1}, ...
                        currTrialNums(iTrial), stimTiming, trialDuration, odorName, concentration, ...
                        flowRate);
            end
        end%iTrial
    end
    
    % Export to .csv file
    odorStims.export_csv(saveDir, 'fileNamePrefix', expMetadata.expID{1})
    
end%if 
catch ME; rethrow(ME); end

%% ----- Flailing epochs -----

moveThresh = 0.05
flowYLim = [0 .2];

try 
figure(1);clf;
set(gcf, 'color', [1 1 1])

for iTrial = 1:expMetadata.nTrials
    
    td = mD(iTrial);  
    
    subaxis(expMetadata.nTrials, 1, iTrial, 'S', 0.02);

    % Optic flow data to identify flailing
    currFlow = meanFlowMags{iTrial};
    currFlow(end) = 0;
    flowFrameDur = median(diff(td.ftFrameTimes));
    flowFrameTimes = (1:1:numel(currFlow)) * flowFrameDur; % NOTE: this is only an  approximation
    plotData = repeat_smooth(currFlow, 20, 'dim', 2, 'smwin', 6);
    plotData = plotData - min(plotData);
    plot(flowFrameTimes, plotData, 'color', 'k');
    hold on;
    if td.usingOptoStim
        hold on; plot_stim_shading([td.optoStimOnsetTimes, td.optoStimOffsetTimes])
    end
    plot([td.ftFrameTimes(1), td.ftFrameTimes(end)], [moveThresh, moveThresh],...
            'linewidth', 0.5, 'color', 'r');
    ylim(flowYLim)
    xlim([0 max(trialMetaData.trialDuration)])
    ylabel('Optic flow (flailing)')
    
end
catch ME; rethrow(ME); end
%%

ftFrameTimes = {ftData.frameTimes};
flowFrameTimes = {};
for iTrial = 1:numel(meanFlowMags)
    
    % NOTE: this is an approximation, but there is mostly a 1:1 correspondence between FicTrac data
    % frames and frames in the captured video, so the median IFI in the FicTrac data should represent
    % the true video frame rate
    flowFrameDur = median(diff(ftFrameTimes{iTrial}));
    flowFrameTimes{iTrial} = (1:1:numel(meanFlowMags{iTrial})) * flowFrameDur; 
    
    % Warn user if there's a discrepency in the flow frame times based on total trial duration
    discrepVal = abs(flowFrameTimes{iTrial}(end) - expMetadata.trialDuration);
    if discrepVal > 0.5
        warning([num2str(discrepVal, 4), ' sec discrepancy in estimated trial duration from flow', ...
            ' frame times for trial #', num2str(iTrial)]);
    end
end

flailingEvents = flailingEvent();
flailingEvents = flailingEvents.append_flow_data(expMetadata.expID{1}, trialMetadata.trialNum, ...
        meanFlowMags, flowFrameTimes, moveThresh);

flailingEvents.export_csv(saveDir, 'fileNamePrefix', expMetadata.expID{1});








