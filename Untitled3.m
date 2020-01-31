
% Identify epochs of quiescence vs. flailing
moveThresh = 0.03;

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

trialNums = [1:10];
panelsStimPositions = 14;%[19 43 67]; % 
quiescenceWinSec = [4 4];

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

roiNum = 2;
smWin = 3;
singleTrials = 1;
skipTimes =[4.02 4.62];


for iPos = 1:numel(stimWinFlData)
    test2 = squeeze(stimWinFlData(iPos).dff(:, roiNum, :));% - squeeze(stimWinFlData(iPos).rawFl(:, 3, :));
%     test2 = repeat_smooth(test2, 10, 'smwin', 3);
    plotTimes = ((1:size(test2, 1)) ./ mD(1).volumeRate) - (1/mD(1).volumeRate);
    
    if ~isempty(skipTimes)
        test2(plotTimes > skipTimes(1) & plotTimes < skipTimes(2), :) = nan;
        plotTimes(plotTimes > skipTimes(1) & plotTimes < skipTimes(2)) = nan;
    end
    
    figure(iPos);clf
    % Plot single trial data
    if singleTrials
        plot(plotTimes, smoothdata(test2, 'gaussian', smWin, 'omitnan'))
    end
    hold on
    
    sd = std(test2, [], 2, 'omitnan');
    avg = mean(smoothdata(test2, 'gaussian', smWin, 'omitnan'), 2); 
    sd = sd(~isnan(plotTimes));
    avg = avg(~isnan(plotTimes));
    sem = sd ./ (size(test2, 2)^0.5);
    jbfill(plotTimes(~isnan(plotTimes)), [avg + sem]', [avg - sem]', 'k', 'k', 1, 0.2);
    
    plot(plotTimes, mean(smoothdata(test2, 'gaussian', smWin, 'omitnan'), 2), 'linewidth', 2, 'color', 'k')
    hold on
    stimDur = stimOffsetTimes(1) - stimOnsetTimes(1);
    plot_stim_shading([quiescenceWinSec(1), quiescenceWinSec(1) + stimDur])
%     xlim([plotTimes(1) + 2, plotTimes(end) - 2])
    hold on;
    

    
    
end








