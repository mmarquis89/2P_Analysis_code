
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\20191119-2_38A11_ChR_60D05_7f\ProcessedData';

% Load ROI data
load(fullfile(parentDir, 'roiData_reg.mat'), 'allROIData');

% Load metadata files
load(fullfile(parentDir, 'daqData.mat'), 'expDaqData');
load(fullfile(parentDir, 'metadata.mat'), 'expMetadata');
load(fullfile(parentDir, 'imagingMetadata.mat'), 'imagingMetadata');

%% Strip out any trials that are not included in the ROI data
expDaqData = expDaqData(ismember([expDaqData.trialNum], [allROIData.trialNum]));
expMetadata = expMetadata(ismember([expMetadata.trialNum], [allROIData.trialNum]));
imagingMetadata = imagingMetadata(ismember([imagingMetadata.trialNum], [allROIData.trialNum]));

% Create struct for keeping track of various calculated parameters
analysisMetadata = []; 
for iTrial = 1:numel(allROIData)
    analysisMetadata(iTrial).trialNum = allROIData(iTrial).trialNum
end

nTrials = numel(allROIData)

%% Calculate volume times in seconds

for iTrial = 1:numel(allROIData)
    analysisMetadata(iTrial).volumeRate = imagingMetadata(iTrial).SI.hRoiManager.scanVolumeRate;
    analysisMetadata(iTrial).nVolumes = imagingMetadata(iTrial).SI.hFastZ.numVolumes;
    analysisMetadata(iTrial).volTimes = (1:1:analysisMetadata(iTrial).nVolumes) / ...
            analysisMetadata(iTrial).volumeRate;
end

%% Calculate all opto stim onset and offset times
for iTrial = 1:nTrials
    aM = analysisMetadata(iTrial);
    aM.trialDuration = expMetadata(iTrial).trialDuration;
    
    currStimTimes = expMetadata(iTrial).optoStimTiming;
    if expMetadata(iTrial).usingOptoStim
        
        
        % Calculate all opto stim onset and offset times
        stimVols = round(currStimTimes * aM.volumeRate);
        
        analysisMetadata(iTrial).stimOnsetTimes = ...
                (currStimTimes(1):sum(currStimTimes(2:3)):aM.trialDuration - currStimTimes(2))';
        analysisMetadata(iTrial).stimOffsetTimes = ...
                (sum(currStimTimes(1:2)):sum(currStimTimes(2:3)):aM.trialDuration)';

%         figure(iTrial);clf; hold on
%         currOptoStimData = expDaqData(iTrial).outputData(:, strcmp(expDaqData(iTrial).columnLabels.out, ...
%                 'OptoStimCommand'));
%         sampTimes = (1:1:numel(currOptoStimData)) / (size(expDaqData(iTrial).outputData, 1) / aM.trialDuration);
%         plot(sampTimes, currOptoStimData);
%         plot_stim_shading([analysisMetadata(iTrial).stimOnsetTimes, analysisMetadata(iTrial).stimOffsetTimes]);
%         title(num2str(allROIData(iTrial).trialNum))
    end
    

    
end

%% Plot entire trial's data for each ROI

currTrial = 1;
smWin = 2;
nReps = 10;

currTrialInd = find([allROIData.trialNum] == currTrial);
currROIData = allROIData(currTrialInd).roiDefs;
aM = analysisMetadata(currTrialInd);

for iROI = 1:numel(currROIData)
    
    figure(iROI); clf; hold on
    plot(aM.volTimes, repeat_smooth(currROIData(iROI).data, nReps, 'Dim', 2, 'SmWin', smWin), ...
        'linewidth', 1)
    plot_stim_shading([aM.stimOnsetTimes, aM.stimOffsetTimes]);
    title(['Trial ', num2str(currTrial), '  ', regexprep(currROIData(iROI).name, '_', '\\_')])
    xlim([0, expMetadata(currTrialInd).trialDuration])
    
end

%% Plot averaged responses to opto stims throughout trial

currTrial = 3;
smWin = 2;
nReps = 15;

% yL = [75 120];

baselineDur = 0.5;
postStimDur = 0.5;

currTrialInd = find([allROIData.trialNum] == currTrial);
currROIData = allROIData(currTrialInd).roiDefs;
aM = analysisMetadata(currTrialInd);
if expMetadata(currTrialInd).usingOptoStim
    
    % Convert plotting window timing from seconds to volumes
    stimDur = expMetadata(currTrialInd).optoStimTiming(2);
    baselineDurVols = floor(baselineDur * aM.volumeRate);
    stimDurVols = floor(stimDur * aM.volumeRate);
    postStimDurVols = round(postStimDur * aM.volumeRate);
    
    stimWinData = [];
    for iROI = 1:numel(currROIData)
        
        % Extract data from periods around opto stims
        for iStim = 1:numel(aM.stimOnsetTimes)
            stimOnsetVol = round(aM.stimOnsetTimes(iStim) * aM.volumeRate);
            stimWinVols = (stimOnsetVol - baselineDurVols):(stimOnsetVol + stimDurVols + ...
                postStimDurVols);
            stimWinData(iROI, iStim, :) = currROIData(iROI).data(stimWinVols);
        end
        
%         stimWinData = smoothdata(stimWinData, 3, 'gaussian', smWin);
        stimWinData = repeat_smooth(stimWinData, nReps, 'Dim', 3, 'SmWin', smWin);
        
        figure(iROI); clf; hold on
        cm = jet(size(stimWinData, 2));
        
        % Plot raw traces
        preStimTimes = (-(baselineDur*aM.volumeRate):1:0) / aM.volumeRate;
        postStimTimes = (1:1:((stimDur + postStimDur)*aM.volumeRate)) / aM.volumeRate;
        plotTimes = [preStimTimes, postStimTimes];
        for iStim = 1:size(stimWinData, 2)
            plot(plotTimes, squeeze(stimWinData(iROI, iStim, :)), 'Color', cm(iStim, :));
        end
        
        % Plot trial-averaged trace
        plot(plotTimes, squeeze(mean(stimWinData(iROI, :, :), 2)), 'Color', 'k', 'linewidth', 2);

        % Add title and stim period shading
        plot_stim_shading([0 stimDur])
        title(['Trial ', num2str(currTrial), '  ', regexprep(currROIData(iROI).name, '_', '\\_')])
%         ylim(yL);
        xlabel('Time from stim onset (s)')
        ylabel('Raw F');
        
        
%         % Save
%         if strcmp(currROIData(iROI).name, 'EB_BU-bi')
%             h = gcf;
%             save_figure(h, fullfile(parentDir, 'figures'), ['trial_', num2str(currTrial), '_EB-BU-bi']);
%         end
%         
        
    end
    
end

%%


avgStimData = squeeze(mean(stimWinData, 2)); % --> [ROI, volume]
stimDataOffset = avgStimData - repmat(min(avgStimData, [], 2), 1, size(avgStimData, 2));
stimDataNorm = stimDataOffset ./ repmat(max(stimDataOffset, [], 2), 1, size(avgStimData, 2));


figure(9);clf; 
imagesc(avgStimData); 
hold on
yL = ylim();
plot([baselineDurVols, baselineDurVols], yL, 'color', 'g', 'linewidth', 3)
plot([baselineDurVols + stimDurVols, baselineDurVols + stimDurVols], yL, 'color', 'r', 'linewidth', 3)
ylim(yL);

figure(10);clf; 
imagesc(stimDataOffset)
hold on
plot([baselineDurVols, baselineDurVols], yL, 'color', 'g', 'linewidth', 3)
plot([baselineDurVols + stimDurVols, baselineDurVols + stimDurVols], yL, 'color', 'r', 'linewidth', 3)
ylim(yL);
ylabel('PB glomerulus')
xlabel('Imaging volumes');
title('Trial 4 - raw F around opto stim');

figure(11);clf; 
imagesc(stimDataNorm)
hold on
plot([baselineDurVols, baselineDurVols], yL, 'color', 'g', 'linewidth', 3)
plot([baselineDurVols + stimDurVols, baselineDurVols + stimDurVols], yL, 'color', 'r', 'linewidth', 3)
ylim(yL);
ylabel('PB glomerulus')
xlabel('Imaging volumes');
title('Trial 4 - normalized raw F around opto stim');

%%

posFunTimes = (1:1:(120*50)) / 50;
posFunRep = repmat(posFun, 1, ceil(120/(numel(posFun)/50)));

posFunTrial = posFunRep(1:numel(posFunTimes));


xPosTimes = (1/(numel(xPos)/120)):(1/(numel(xPos)/120)):trialDuration;

xPos = xPos - min(xPos);

figure(1);clf;hold on;
plot(posFunTimes, posFunTrial);
plot(xPosTimes, xPos * 10)











