
% Identify epochs of quiescence vs. flailing
moveThresh = 0.08;

for iTrial = 1:numel(mD)

   currTrialFlow = meanFlowMags{iTrial};
   smFlow = repeat_smooth(currTrialFlow, 20, 'dim', 2, 'smwin', 6);
   smFlow = smFlow - min(smFlow(:));
   moveFrames = smFlow > moveThresh;
   
   % For each quiescence frame, find the distance to the nearest movement frame
   moveFrameDist = zeros(1, numel(moveFrames));
   moveFrameInds = find(moveFrames);
   for iFrame = 1:numel(moveFrames)
       moveFrameDist(iFrame) = min(abs(moveFrameInds - iFrame));
   end
   
   mD(iTrial).moveFlowFrames = moveFrames;
   mD(iTrial).moveFlowFrameDist = moveFrameDist;
   
end



%% Extract fluorescence data for visual stim epochs only when fly is quiescent

trialNums = [5:13];
panelsStimPositions = 14;% [19 67 43]; %
stimPosNames = {'full field flash'}; % {'left bar', 'right bar', 'center bar'}; %
quiescenceWinSec = [6 6];

stimWinFlData = struct();
for iTrial = 1:numel(trialNums)
    
    currTrialInd = find([mD.trialNum] == trialNums(iTrial));
    
    % Identify stim epochs for each panels stim position
    currBarPosData = mD(currTrialInd).panelsPosX;
    stimOnsetFrames = {}; stimOffsetFrames = {};
    for iPos = 1:numel(panelsStimPositions)
        stimOnsetFrames = [];
        stimOffsetFrames = [];
        if iTrial == 1
            stimWinFlData(iPos).rawFl = [];
            stimWinFlData(iPos).dff = [];
            stimWinFlData(iPos).expDff = [];
            stimWinFlData(iPos).zscore = [];
        end
        for iFrame = 2:(numel(currBarPosData)-1)
           if currBarPosData(iFrame - 1) ~= panelsStimPositions(iPos) && ...
                        currBarPosData(iFrame) == panelsStimPositions(iPos)
              stimOnsetFrames(end + 1) = iFrame;
           elseif currBarPosData(iFrame + 1) ~= panelsStimPositions(iPos) && ...
                        currBarPosData(iFrame) == panelsStimPositions(iPos)
              stimOffsetFrames(end + 1) = iFrame;      
           end
        end
        
        % Convert from panels frames to imaging volumes
        stimOnsetTimes = mD(currTrialInd).panelsFrameTimes(stimOnsetFrames);
        stimOffsetTimes = mD(currTrialInd).panelsFrameTimes(stimOffsetFrames);
        
        
        % Filter out stim epochs that are accompanied by movement
        currFlow = meanFlowMags{currTrialInd};
        currFlow(end) = 0;
        flowFrameDur = median(diff(mD(currTrialInd).ftFrameTimes));
        flowFrameTimes = (1:1:numel(currFlow)) * flowFrameDur; % NOTE: this is only an approximation
        
        stimOnsetFlowFrames = []; stimOffsetFlowFrames = [];
        for iStim = 1:numel(stimOnsetTimes)
           [~, stimOnsetFlowFrames(iStim)] = min(abs(stimOnsetTimes(iStim) - flowFrameTimes));
           [~, stimOffsetFlowFrames(iStim)] = min(abs(stimOffsetTimes(iStim) - flowFrameTimes));
        end
        quiescenceWinFlowFrames = ceil(quiescenceWinSec ./ flowFrameDur);
        stimWinStartFlowFrames = stimOnsetFlowFrames - quiescenceWinFlowFrames(1);
        stimWinEndFlowFrames = stimOffsetFlowFrames + quiescenceWinFlowFrames(2);
        
        stimOnsetTimes(stimWinStartFlowFrames < 1) = [];
        stimOffsetTimes(stimWinStartFlowFrames < 1) = [];
        stimWinEndFlowFrames(stimWinStartFlowFrames < 1) = [];
        stimWinStartFlowFrames(stimWinStartFlowFrames < 1) = [];
        
        stimOnsetTimes(stimWinEndFlowFrames > ...
                numel(mD(currTrialInd).moveFlowFrameDist)) = [];
        stimOffsetTimes(stimWinEndFlowFrames > ...
                numel(mD(currTrialInd).moveFlowFrameDist)) = [];
        stimWinEndFlowFrames(stimWinEndFlowFrames > ...
                numel(mD(currTrialInd).moveFlowFrameDist)) = [];
        stimWinStartFlowFrames(stimWinEndFlowFrames > ...
                numel(mD(currTrialInd).moveFlowFrameDist)) = [];
            
        minMoveDists = [];
        for iStim = 1:numel(stimOnsetTimes)
            minMoveDists(iStim) = min(mD(currTrialInd).moveFlowFrameDist( ...
                    stimWinStartFlowFrames(iStim):stimWinEndFlowFrames(iStim)));
        end
        
        stimWinOnsetTimes = stimOnsetTimes(logical(minMoveDists)) - quiescenceWinSec(1);
        stimWinOffsetTimes = stimOffsetTimes(logical(minMoveDists)) + quiescenceWinSec(2);
         
        % Extract fluorescence data from target stim window for each ROI
        stimWinOnsetVols = floor(stimWinOnsetTimes .* mD(currTrialInd).volumeRate);
        stimWinOffsetVols = floor(stimWinOffsetTimes .* mD(currTrialInd).volumeRate);
        for iStim = 1:numel(stimWinOnsetVols)
            currTd = mD(currTrialInd);
            currStart = stimWinOnsetVols(iStim);
            currEnd = stimWinOffsetVols(iStim);
            if ~isempty(stimWinFlData(iPos).rawFl) && ...
                        (currEnd - currStart + 1) > size(stimWinFlData(iPos).rawFl, 1)
               currEnd = currEnd - 1;
            elseif ~isempty(stimWinFlData(iPos).rawFl) && ...
                        (currEnd - currStart + 1) < size(stimWinFlData(iPos).rawFl, 1)
                currEnd = currEnd + 1;
            end
            stimWinFlData(iPos).rawFl = cat(3, stimWinFlData(iPos).rawFl, ...
                    currTd.rawFlMat(currStart:currEnd, :));     % --> [vol, ROI, stim]
            stimWinFlData(iPos).dff = cat(3, stimWinFlData(iPos).dff, ...
                    currTd.dffMat(currStart:currEnd, :));       % --> [vol, ROI, stim]
            stimWinFlData(iPos).expDff = cat(3, stimWinFlData(iPos).expDff, ...
                    currTd.expDffMat(currStart:currEnd, :));    % --> [vol, ROI, stim]
            stimWinFlData(iPos).zscore = cat(3, stimWinFlData(iPos).zscore, ...
                    currTd.zscoreMat(currStart:currEnd, :));    % --> [vol, ROI, stim]
        end
        
        
    end%iPos    
    
    
end%iTrial

for iPos = 1:numel(stimWinFlData)
   disp(['Stim position ', num2str(iPos), ': ', num2str(size(stimWinFlData(iPos).rawFl, 3)), ...
            ' presentations without fly movement']); 
end
disp(' ')

%%

saveFig = 0;

roiNum = 1;
smWin = 1;
singleTrials = 1;
skipTimes =[-0.05 0.45];

dataType = 'dff';


% Adjust figure dimensions for number of subplots
if numel(stimWinFlData) == 1
   figDims = [800, 500];
   mr = 0.05;
   ml = 0.12;
else
   figDims = [numel(stimWinFlData) * 600, 500];
   mr = 0.03;
   ml = 0.06;
end

% Create figure
f = figure(1);clf
f.Color = [1 1 1];
f.Position(3:4) = figDims;

roiName = mD(trialNums(1)).roiData(roiNum).name;
clear ax
for iPos = 1:numel(stimWinFlData)
    plotData = squeeze(stimWinFlData(iPos).(dataType)(:, roiNum, :));% - squeeze(stimWinFlData(iPos).rawFl(:, 3, :));
    plotTimes = ((1:size(plotData, 1)) ./ mD(1).volumeRate) - (1/mD(1).volumeRate) ...
            - quiescenceWinSec(1);

    if ~isempty(skipTimes)
        plotData(plotTimes > skipTimes(1) & plotTimes < skipTimes(2), :) = nan;
        plotTimes(plotTimes > skipTimes(1) & plotTimes < skipTimes(2)) = nan;
    end

    ax(iPos) = subaxis(1, numel(stimWinFlData), iPos, ...
            'mb', 0.14, 'mt', 0.15, 'mr', mr, 'ml', ml);
    hold on
    
    % Plot single trial data
    cm = jet(size(plotData, 2)) * 0.95;
    if singleTrials
        for iTrial = 1:size(plotData, 2)
            plot(plotTimes, smoothdata(plotData(:, iTrial), 'gaussian', smWin, 'omitnan'), ...
                    'color', cm(iTrial, :), 'linewidth', 1)
        end
    end
    
    % Shade +/- SEM
    sd = std(plotData, [], 2, 'omitnan');
    avg = mean(smoothdata(plotData, 'gaussian', smWin, 'omitnan'), 2); 
    sd = sd(~isnan(plotTimes));
    avg = avg(~isnan(plotTimes));
    sem = sd ./ (size(plotData, 2)^0.5);
    jbfill(plotTimes(~isnan(plotTimes)), [avg + sem]', [avg - sem]', 'k', 'k', 1, 0.2);
    
    % Plot mean
    plot(plotTimes, mean(smoothdata(plotData, 'gaussian', smWin, 'omitnan'), 2), ...
            'linewidth', 3, 'color', 'k')
    hold on
    
    % Shade stim period
    stimDur = stimOffsetTimes(1) - stimOnsetTimes(1);
    plot_stim_shading([0 stimDur])
    hold on;
    
    % Add labels to plot
    xlim([plotTimes(1), plotTimes(end)]);
    xlabel('Time from stim onset (s)');
    ax(iPos).FontSize = 14;
    if strcmp(dataType, 'expDff')
        yLabStr = 'Full exp. dF/F';
    elseif strcmp(dataType, 'dff')
        yLabStr = 'Trial dF/F';
    elseif strcmp(dataType, 'zscore')
        yLabStr = 'z-scored dF/F';
    elseif strcmp(dataType, 'rawFl')
        yLabStr = 'Raw Fl';
    else
        ylabel('');
    end
    if iPos == 1
       ylabel(yLabStr); 
    end
    
%     title({[mD(1).expID, '  -  mean ', roiName, ' ', stimPosNames{iPos}, ' responses'], ...
%         ' (no fly movement within plot epochs)'}, 'fontsize', 14)
    title(stimPosNames{iPos}, 'fontsize', 14)
end

linkaxes(ax);

suptitle([[mD(1).expID, '  -  mean ', roiName, ' visual stim responses'], ...
    ' (no fly movement within plot epochs)']);

%%

winSize = 2;



stimStartVol = find(isnan(plotTimes), 1, 'first');
stimEndVol = find(isnan(plotTimes), 1, 'last');
winSizeVols = round(winSize * mD(1).volumeRate);

preStimVols = (stimStartVol - winSizeVols):(stimStartVol - 1);
postStimVols = (stimEndVol + 1):(stimEndVol + winSizeVols);

preStimMean = mean(plotData(preStimVols, :), 1);
postStimMean = mean(plotData(postStimVols, :), 1);

f = figure(2);clf; 
f.Color = [1 1 1];
hold on

plot([ones(1, numel(preStimMean)); 2 * ones(1, numel(preStimMean))], ...
         [preStimMean; postStimMean], '-o', ...
         'linewidth', 1);
xlim([0.5 2.5])

plot([mean(preStimMean), mean(postStimMean)], '-s', 'linewidth', 3, ...
        'markersize', 10, 'color', 'k', 'markerfacecolor', 'k')

[~, p] = ttest(preStimMean, postStimMean)





