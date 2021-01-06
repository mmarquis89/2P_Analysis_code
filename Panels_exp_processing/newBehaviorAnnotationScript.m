
startDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';

expDir = uigetdir(startDir, 'Select an experiment directory');
outputDir = fullfile(expDir, 'ProcessedData');

if ~isfolder(outputDir)
    mkdir(outputDir)
end

expID = regexp(expDir, '(?<=\\)........-.(?=_)', 'match');

%%

flyFlowThresh = 0.05;
ballFlowThresh = 0.03;

smWin = 5;
smReps = 6;

minIsoMoveDur = 0.25;
minLocEpochDur = 0.25;
minQuiescenceDur = 0.25;

try 
    
% Load data
ftData = load_ft_data(expID, outputDir);
[expMd, trialMd] = load_metadata(expID, outputDir);

% Get flow data
vidFrameTimes = ftData.vidFrameTimes;
flyFlow = ftData.meanFlow;
ballFlow = ftData.meanFlowBall;
nTrials = numel(vidFrameTimes);

% Preprocess data and generate state maps
globalMins = [inf, inf];
globalMaxes = [-inf, -inf];
for iTrial = 1:nTrials
    smFlyFlow = repeat_smooth(flyFlow{iTrial}, smReps, 'smWin', smWin);
    smBallFlow = repeat_smooth(ballFlow{iTrial}, smReps, 'smWin', smWin);
    globalMins(1) = min([globalMins(1), min(smFlyFlow)]);
    globalMins(2) = min([globalMins(2), min(smBallFlow)]);
    globalMaxes(1) = max([globalMaxes(1), max(smFlyFlow)]);
    globalMaxes(2) = max([globalMaxes(2), max(smBallFlow)]);
end

flyFlowNorm = {};
ballFlowNorm = {};
combThreshData = {};
for iTrial = 1:nTrials
    
    % Smooth data
    smFlyFlow = repeat_smooth(flyFlow{iTrial}, smReps, 'smWin', smWin);
    smBallFlow = repeat_smooth(ballFlow{iTrial}, smReps, 'smWin', smWin);
    
    % Scale from 0-1
    quick_norm = @(x, minVal, maxVal) (x - minVal)  ./  (maxVal - minVal);    
    flyFlowNorm{iTrial} = quick_norm(smFlyFlow, globalMins(1), globalMaxes(1));
    ballFlowNorm{iTrial} = quick_norm(smBallFlow, globalMins(2), globalMaxes(2));
    
    % Calculate thresholded state maps
    thresholdedFlyFlow = flyFlowNorm{iTrial} > flyFlowThresh;
    thresholdedBallFlow = ballFlowNorm{iTrial} > ballFlowThresh;
    combThreshData{iTrial} = thresholdedFlyFlow + (2 * thresholdedBallFlow);
    
    % Set the first and last frames to quiescence for all trials
    combThreshData{iTrial}(1) = 0;
    combThreshData{iTrial}(end) = 0;
    
end


% --------- CLEANUP --------------------------------------------------------------------------------

newCombThresh = combThreshData;

for iTrial = 1:nTrials
    
    currCombThresh = newCombThresh{iTrial};
    nFrames = numel(currCombThresh);
    currFrameTimes = vidFrameTimes{iTrial};
    
    % Convert Loc+, Move- epochs to quiescence
    currCombThresh(currCombThresh == 2) = 0;
    
    % Eliminate any too-short locomotion bouts
    currInd = 1;
    currVal = currCombThresh(currInd);
    while currInd < nFrames
        
        % Move until currVal is Loc+
        while currVal ~= 3 && currInd < nFrames
            lastVal = currVal;
            locStartInd = currInd;
            currInd = currInd + 1;
            currVal = currCombThresh(currInd);
        end
        
        % Move to end of epoch and count length
        while currVal == 3 && currInd < nFrames
            currInd = currInd + 1;
            currVal = currCombThresh(currInd);
        end
        
        % If epoch is too short convert it to preceding action
        if (currFrameTimes(currInd) - currFrameTimes(locStartInd)) < minLocEpochDur
            currCombThresh(locStartInd:currInd) = lastVal;
        end
    end
    
    % Trim Move+, Loc- epochs from the start and end of locomotion bouts
    currInd = 1;
    postLoc = 0;
    currVal = currCombThresh(currInd);
    while currInd < nFrames
        
        % Move until currVal is Move+
        while currVal ~= 1 && currInd < nFrames
            if currVal == 3
                postLoc = 1;
            else
                postLoc = 0;
            end
            currInd = currInd + 1;
            currVal = currCombThresh(currInd);
        end
        
        if currInd < nFrames
            
            % Mark isoMove start frame
            isoMoveStartInd = currInd;
            
            % Move to end of isoMove epoch
            while currVal == 1 && currInd < nFrames
                currInd = currInd + 1;
                currVal = currCombThresh(currInd);
            end
            epochDur = currFrameTimes(currInd - 1) - currFrameTimes(isoMoveStartInd);
            
            % Remove Loc-, Move+ epoch if necessary
            if epochDur < minIsoMoveDur
                if postLoc
                    if currVal == 0
                        % Remove trailing Loc-, Move+ epoch
                        currCombThresh(isoMoveStartInd:currInd) = 0;
                        postLoc = 0;
                    elseif currVal == 3
                        % Convert sandwich Loc-, Move+ epoch to Loc+, Move+
                        currCombThresh(isoMoveStartInd:currInd) = 3;
                    end
                else
                    if currVal == 3
                        % Remove preceding Loc-, Move+ epoch
                        currCombThresh(isoMoveStartInd:currInd) = 0;
                    end
                end
            end
            
        end%if
    end%while
    
    
    % Eliminate any too-short quiescence sandwiches
    currInd = 1;
    currVal = currCombThresh(currInd);
    while currInd < nFrames
        
        % Move until currVal is Loc-, Move-
        while currVal ~= 0 && currInd < nFrames
            lastVal = currVal;
            locStartInd = currInd;
            currInd = currInd + 1;
            currVal = currCombThresh(currInd);
        end
        
        % Move to end of epoch and count length
        while currVal == 0 && currInd < nFrames
            currInd = currInd + 1;
            currVal = currCombThresh(currInd);
        end
        
        % If epoch is too short convert it to preceding action
        if (currFrameTimes(currInd) - currFrameTimes(locStartInd)) < minQuiescenceDur
            currCombThresh(locStartInd:currInd) = lastVal;
        end
    end
    
    newCombThresh{iTrial} = currCombThresh;

end%iTrial

catch ME; rethrow(ME); end

% --------- VISUALIZE PRE- AND POST- CLEANUP -------------------------------------------------------

trialNum = 1;

% Imagesc before and after comparison
f = figure(1); clf;
f.Color = [1 1 1];
clear ax;
ax(1) = subaxis(4, 1, 1, 'ml', 0.02, 'mr', 0.02, 'sv', 0.05, 'mb', 0.08);
imagesc([vidFrameTimes{trialNum}(1), vidFrameTimes{trialNum}(end)], [0, 1], [combThreshData{trialNum}']);
colormap(gca, [rgb('Indigo'); rgb('Orange'); rgb('Yellow'); rgb('Cyan')]);
title('Thresholded')
ax(1).YTickLabel = [];
ax(1).XTickLabel = [];
ax(1).FontSize = 14;

ax(2) = subaxis(4, 1, 2); hold on;
imagesc([vidFrameTimes{trialNum}(1), vidFrameTimes{trialNum}(end)], [0, 1], [newCombThresh{trialNum}']);
colormap(gca, [rgb('Indigo'); rgb('Orange'); rgb('Yellow'); rgb('Cyan')]);
title('Final')
ax(2).YTickLabel = [];
ax(2).FontSize = 14;
linkaxes(ax(1:2), 'x');
xlim([vidFrameTimes{trialNum}(1), vidFrameTimes{trialNum}(end)])

% Overlay before and after comparison
ax(3) = subaxis(4, 1, [3 4]); hold on;
xx = [vidFrameTimes{trialNum}'];
plot(xx, [flyFlowNorm{trialNum}'], 'linewidth', 1);
plot(xx, [ballFlowNorm{trialNum}'], 'linewidth', 1);
plot(xx, [newCombThresh{trialNum}'] ./ 3, 'linewidth', 2);
plot([xx(1), xx(end)], [1, 1] * flyFlowThresh, 'color', 'r', 'linewidth', 2);
plot([xx(1), xx(end)], [1, 1] * ballFlowThresh, '--', 'color', 'g', 'linewidth', 2);
legend('Fly flow', 'Ball flow', 'FinalAnnot', 'Fly flow thresh', 'Ball flow thresh', ...
        'autoupdate', 'off');
ylim([0 1]);
xlim([xx(1), xx(end)]);
ax(3).FontSize = 14;


%% Create and save behavior event objects for locomotion and isolated movement

isoMoveEvents = behaviorEvent('IsolatedMovement');
locEvents = behaviorEvent('Locomotion');
for iTrial = 1:nTrials
    currAnnotData = newCombThresh{iTrial};
    currFrameTimes = vidFrameTimes{iTrial};
    
    % Separate according to behavior type
    iFrames = currAnnotData == 1;
    lFrames = currAnnotData == 3;
    
    if sum(iFrames) > 0
        isoMoveEvents = isoMoveEvents.append_annotation_data(expID{:}, ftData.trialNum(iTrial), ...
                iFrames, currFrameTimes);
    end
    if sum(lFrames) > 0
        locEvents = locEvents.append_annotation_data(expID{:}, ftData.trialNum(iTrial), ...
                lFrames, currFrameTimes);
    end
end

if ~isempty(isoMoveEvents.eventData)
    isoMoveEvents.export_csv(outputDir, 'fileNamePrefix', expID{:});
end
if ~isempty(locEvents.eventData)
    locEvents.export_csv(outputDir, 'fileNamePrefix', expID{:});
end

p = [];
p.flyFlowThresh = flyFlowThresh;
p.ballFlowThresh = ballFlowThresh;
p.smWin = smWin;
p.smReps = smReps;
p.minIsoMoveDur = minIsoMoveDur;
p.minLocEpochDur = minLocEpochDur;
p.minQuiescenceDur = minQuiescenceDur;

save(fullfile(outputDir, 'behaviorAnnotationParams.mat'), 'p')

disp('Behavior annotations saved');