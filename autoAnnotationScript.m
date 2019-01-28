% ==================================================================================================
% LOAD DATA
% ==================================================================================================

expDate = '2019_01_24_exp_1';
sid = 0;
FRAME_RATE = 25;
trialDuration = 20;

behaviorLabels = {'Quiescence', 'Locomotion', 'IsolatedMovement'};

imgDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)'];
vidDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
imgDir = vidDir;

% Get frame count info
frameInfo = [];
[goodTrials, frameCounts, allTrialsFrameCount, nFrames] = frame_count_check(vidDir, sid, FRAME_RATE, ...
    trialDuration);
frameTimes = (1:(FRAME_RATE*trialDuration)) * (1/FRAME_RATE); 
frameInfo.goodTrials = goodTrials;
frameInfo.nFrames = nFrames;
frameInfo.frameTimes = frameTimes;
frameInfo.FRAME_RATE = FRAME_RATE;
nTrials = numel(goodTrials);


% % Make sure frame counts are consistent with each other
% assert(sum(frameCounts) == allTrialsFrameCount, ...
%     'Error: sum of individual frame counts is not equal to concatenated video frame count');

% Load FicTrac data
ftDir = fullfile(vidDir, 'FicTracData');
ftData = load_fictrac_data(frameInfo, 'Sid', sid, 'ParentDir', ftDir);
goodTrials(ftData.badFtTrials) = 0;

% Load flow data
load(fullfile(vidDir, ['sid_', num2str(sid), '_flow_data_norm.mat'])) % flyFlowNorm
flowArr = [];
for iTrial = 1:nTrials
    if goodTrials(iTrial)
        flowArr(:, iTrial) = flyFlowNorm{iTrial}'; % --> [frame, trial]
    else
        flowArr(:, iTrial) = nan(1, nFrames); % Fill in bad trials with nan placeholders
    end
end
flowArr(1,:) = flowArr(2, :); % Ignore artificially high first frame of each trial


% Get old-format annot data
try
load(fullfile(imgDir, 'Annotations.mat'), 'trialAnnotations'); % trialAnnotations
annotArr = [];
for iTrial = 1:nTrials
    if goodTrials(iTrial)
        annotArr(:,iTrial) = trialAnnotations{iTrial}.actionNums;    % --> [frame, trial]
    elseif iTrial == nTrials
        annotArr(:,end + 1) = zeros(size(annotArr, 1), 1);
    end
end
catch
    annotArr = zeros(size(flowArr));
end

%% Save annotation array and annot params

% Save parameters
annotParams = struct('expDate', expDate, 'flowThresh', flowThresh, 'moveThresh', moveThresh, ...
                            'smWin', smWin, 'moveSmReps', moveSmReps, 'flowSmReps', flowSmReps, ...
                            'minIsoMoveLen', minIsoMoveLen, 'minLocEpochLen', minLocEpochLen, ...
                            'minQuiescenceLen', minQuiescenceLen);
                        
if ~isdir(imgDir)
    mkdir(imgDir)
end
                        
save(fullfile(imgDir, 'autoAnnotations.mat'), 'trialAnnotations', 'annotParams', 'ftData', 'flowArr', 'goodTrials', 'behaviorLabels', 'frameInfo')
close all

%%
t = 171;
flowThresh = 0.04;
moveThresh = 0.035;

smWin = 3;
smWinAlt = 1;

moveSmReps = 6;
flowSmReps = 4;

minIsoMoveLen = 6;
minLocEpochLen = 6;
minQuiescenceLen = 6;

% Extract variables from ftData
moveSpeed = ftData.moveSpeed * FRAME_RATE * 4.5;      % --> [frame, trial]
fwSpeed = ftData.fwSpeed * FRAME_RATE * 4.5;          % --> [frame, trial]
yawSpeed = rad2deg(ftData.yawSpeed * FRAME_RATE);     % --> [frame, trial]

% Smooth data
smMoveSpd = repeat_smooth(moveSpeed, smWin, 1, moveSmReps);
smYawSpd = repeat_smooth(yawSpeed, smWin, 1, moveSmReps);
smFlow = repeat_smooth(flowArr, smWin, 1, flowSmReps);
smFlowAlt = repeat_smooth(flowArr, smWinAlt, 1, 1);
smMoveAlt = repeat_smooth(moveSpeed, smWinAlt, 1, 1);

% Scale from 0-1
quick_norm = @(x) (x - min(abs(x(:))))  ./  max( abs(x(:)) - min(abs(x(:))) );
moveNorm = quick_norm(smMoveSpd);
yawNorm = quick_norm(smYawSpd);
flowNorm = quick_norm(smFlow);
flowNormAlt = quick_norm(smFlowAlt);
annotNorm = quick_norm(annotArr);

moveYawAvg = mean(cat(3, moveNorm, abs(yawNorm)), 3);
moveYawMax = max(cat(3, moveNorm, abs(yawNorm)), [], 3);

% Calculate thresholded state maps
thresholdedFlow = flowNorm > flowThresh;
thresholdedMove = moveYawMax > moveThresh;
combThreshData = thresholdedFlow + (2 * thresholdedMove);
annotPlot = cat(1, ones(1,nTrials), 0.75 * ones(1,nTrials), 0.5 * ones(1,nTrials), zeros(1,nTrials), annotNorm);
combThreshPlot = cat(1, zeros(1,nTrials), ones(1,nTrials), 2 * ones(1,nTrials), 3 * ones(1,nTrials), combThreshData);



% % Imagesc plots of entire experiment
% figure(8);clf
% subaxis(2, 1, 1)
% imagesc(annotPlot'); colormap(gca, [rgb('Indigo'); rgb('Magenta'); rgb('Cyan'); rgb('Gold'); rgb('Gold')])
% subaxis(2, 1, 2);
% imagesc(combThreshPlot'); colormap([rgb('Indigo');rgb('Orange');rgb('yellow');rgb('Cyan')]);


% --------- CLEANUP --------------------------------------------------------------------------------

% Set first and last frames to quiescence for all trials
newCombThresh = combThreshData;
newCombThresh(1,:) = 0;
newCombThresh(end,:) = 0;

% Convert Move+, Flow- epochs to quiescence
newCombThresh(newCombThresh == 2) = 0;

% Eliminate any too-short locomotion bouts
for iTrial = 1:nTrials
    currInd = 1;
    currVal = newCombThresh(currInd, iTrial);
    while currInd < nFrames
        
        % Move until currVal is Loc
        while currVal ~= 3 && currInd < nFrames
            lastVal = currVal;
            locStartInd = currInd;
            currInd = currInd + 1;
            currVal = newCombThresh(currInd, iTrial);
        end
        
        % Move to end of epoch and count length
        while currVal == 3 && currInd < nFrames
            currInd = currInd + 1;
            currVal = newCombThresh(currInd, iTrial);
        end
        
        % If epoch is too short convert it to preceding action
        if (currInd - locStartInd) < minLocEpochLen
            newCombThresh(locStartInd:currInd, iTrial) = lastVal;
        end
    end
end

% Trim isoMove epochs from the start and end of locomotion bouts
for iTrial = 1:nTrials
    currInd = 1;
    postLoc = 0;
    currVal = newCombThresh(currInd, iTrial);
    while currInd < nFrames
        
        % Move until currVal is isoMove
        while currVal ~= 1 && currInd < nFrames
            if currVal == 3
                postLoc = 1;
            else
                postLoc = 0;
            end
            currInd = currInd + 1;
            currVal = newCombThresh(currInd, iTrial);
        end
        
        if currInd < nFrames
            
            % Mark isoMove start frame
            isoMoveStartInd = currInd;
            
            % Move to end of isoMove epoch
            while currVal == 1 && currInd < nFrames
                currInd = currInd + 1;
                currVal = newCombThresh(currInd, iTrial);
            end
            epochLen = (currInd - 1) - isoMoveStartInd;
            
            % Remove isoMove epoch if necessary
            if epochLen < minIsoMoveLen
                if postLoc
                    if currVal == 0
                        % Remove trailing isoMove
                        newCombThresh(isoMoveStartInd:currInd, iTrial) = 0;
                        postLoc = 0;
                    elseif currVal == 3
                        % Convert sandwich isoMove to locomotion
                        newCombThresh(isoMoveStartInd:currInd, iTrial) = 3;
                    end
                else
                    if currVal == 3
                        % Remove preceding isoMove
                        newCombThresh(isoMoveStartInd:currInd, iTrial) = 0;
                    end
                end
            end
            
        end%if
    end%while
end%for

% Eliminate any too-short quiescence sandwiches
for iTrial = 1:nTrials
    currInd = 1;
    currVal = newCombThresh(currInd, iTrial);
    while currInd < nFrames
        
        % Move until currVal is quiescence
        while currVal ~= 0 && currInd < nFrames
            lastVal = currVal;
            locStartInd = currInd;
            currInd = currInd + 1;
            currVal = newCombThresh(currInd, iTrial);
        end
        
        % Move to end of epoch and count length
        while currVal == 0 && currInd < nFrames
            currInd = currInd + 1;
            currVal = newCombThresh(currInd, iTrial);
        end
        
        % If epoch is too short convert it to preceding action
        if (currInd - locStartInd) < minQuiescenceLen
            newCombThresh(locStartInd:currInd, iTrial) = lastVal;
        end
    end
end

% Save final array to different variable
trialAnnotations = newCombThresh'; % --> [trial, frame]

% Visualize pre- and post-cleanup
f = figure(9); clf
xt = 0:FRAME_RATE:nFrames;
xtl = 0:trialDuration;
f.Color = [1 1 1]; f.Position = [1921 -623 1080 1843];
newCombThreshPlot = cat(1, zeros(1,nTrials), ones(1,nTrials), 2 * ones(1,nTrials), 3 * ones(1,nTrials), newCombThresh);
subaxis(2, 1, 1)
imagesc(combThreshPlot'); colormap(gca, [rgb('Indigo');rgb('Orange');rgb('yellow');rgb('Cyan')]); title('Thresholded')
% imagesc(annotPlot'); colormap(gca, [rgb('Indigo'); rgb('Magenta'); rgb('Cyan'); rgb('Gold'); rgb('Gold')]) title('Manual')
set(gca, 'XTick', xt); set(gca, 'XTickLabel', xtl);
subaxis(2, 1, 2);
imagesc(newCombThreshPlot'); colormap(gca, [rgb('Indigo');rgb('Orange');rgb('yellow');rgb('Cyan')]); title('Final')
set(gca, 'XTick', xt); set(gca, 'XTickLabel', xtl);

% Overlay
f = figure(5);clf;hold on;
f.Color = [1 1 1]; f.Position = [-1917 48 1916 420];
legendStr = [];
xData = 0:(1/FRAME_RATE):trialDuration-(1/FRAME_RATE);
% plot(xData, moveNorm(:,t)); legendStr{end + 1} = 'Move';
% plot(xData, abs(fwNorm)); legendStr{end + 1} = 'FW';
% plot(xData, abs(yawNorm(:,t))); legendStr{end + 1} = 'Yaw';
plot(xData, flowNorm(:,t)); legendStr{end + 1} = 'Flow';
plot(xData, flowNormAlt(:,t)); legendStr{end + 1} = 'SingleSmFlow';
% plot(xData, moveNormAlt(:,t)); legendStr{end + 1} = 'SingleSmMove';
% plot(xData, annotNorm(:,t), 'linewidth', 2); legendStr{end + 1} = 'Annot';
% plot(xData, moveYawAvg(:,t)); legendStr{end + 1} = 'MoveYawAvg';
plot(xData, moveYawMax(:,t)); legendStr{end + 1} = 'MoveYawMax';
% plot(xData, combThreshPlot(5:end, t)./3, 'linewidth', 2); legendStr{end + 1} = 'AutoAnnot';
plot(xData, newCombThreshPlot(5:end, t)./3, 'linewidth', 2); legendStr{end + 1} = 'CleanAutoAnnot';
legend(legendStr)
ylim([0 1]);

% Imagesc view of old and new annotations
f = figure(7);clf;
f.Color = [1 1 1]; f.Position = [-1863 545 1855 420];
subaxis(3, 1, 1)
imagesc(annotPlot(:,t)'); 
colormap(gca, [rgb('Indigo'); rgb('Magenta'); rgb('Cyan'); rgb('DarkRed'); rgb('Gold')]);
ax = gca(); ax.XTick = []; ax.YTick = []; ylabel(ax, 'Manual');
subaxis(3, 1, 2)
imagesc(combThreshPlot(:,t)');
colormap(gca, [rgb('Indigo');rgb('Orange');rgb('yellow');rgb('Cyan')])
ax = gca(); ax.XTick = []; ax.YTick = []; ylabel(ax, 'Thresholded');
subaxis(3, 1, 3)
imagesc(newCombThreshPlot(:,t)');
colormap(gca, [rgb('Indigo');rgb('Orange');rgb('yellow');rgb('Cyan')]);
ax = gca();
ax.XTick = 0:FRAME_RATE:nFrames;
ax.XTickLabel = 0:trialDuration;
ax.YTick = [];
ylabel(ax, 'Final');








