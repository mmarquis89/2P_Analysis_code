%% COMBINE DATA FROM A BLOCK OF COMPATIBLE TRIALS

blTrials = [1:17];

bD = mD(ismember([mD.trialNum], blTrials));

% Check that have compatible values in scalar fields 
if numel(unique([bD.daqSampRate])) > 1 || ...
        numel(unique({bD.expID})) > 1 || ...
        numel(unique([bD.nDaqSamples])) > 1 || ...
        numel(unique([bD.nPlanes])) > 1 || ...
        numel(unique([bD.nPanelsFrames])) > 1 || ...
        numel(unique([bD.nVolumes])) > 1 || ...
        numel(unique([bD.panelsDisplayRate])) > 1 || ...
        numel(unique([bD.trialDuration])) > 1 || ...
        numel(unique([bD.using2P])) > 1 || ...
        numel(unique([bD.volumeRate])) > 1 || ...
        numel(unique([bD.nVolumes])) > 1 || ...
        numel(unique([bD.panelsCycleFrames])) > 1 || ...
        numel(unique([bD.panelsCycleTime])) > 1
    disp('Error: the specified trials are not compatible');
end

% Add data that applies to all trials
bl = [];
bl.daqSampRate = bD(1).daqSampRate;
bl.daqSampTimes = bD(1).daqSampTimes;
bl.expID = bD(1).expID;
bl.nDaqSamples = bD(1).nDaqSamples;
bl.nPanelsFrames = bD(1).nPanelsFrames;
bl.nPlanes = bD(1).nPlanes;
bl.nVolumes = bD(1).nVolumes;
bl.optoStimTiming = bD(1).optoStimTiming;
bl.panelsCycleFrames = bD(1).panelsCycleFrames;
bl.panelsCycleTime = bD(1).panelsCycleTime;
bl.panelsDisplayRate = bD(1).panelsDisplayRate;
bl.panelsFrameTimes = bD(1).panelsFrameTimes;
bl.panelsPattern = bD(1).panelsPattern;
bl.panelsPosX = bD(1).panelsPosX;
bl.panelsPosY = bD(1).panelsPosY;
bl.refImages = bD(1).refImages;
bl.trialDuration = bD(1).trialDuration;
bl.volTimes = bD(1).volTimes;


% Add compiled data
bl.dffArr = []; bl.dffVectAvgRad = []; bl.dffVectAvgWedge = []; bl.dffVectStrength = [];
bl.ftData = []; bl.optoStimOnsetTimes = []; bl.optoStimOffsetTimes = [];
bl.wedgeDffArr = []; bl.wedgeZscoreArr = []; bl.zscoreArr = []; bl.usingOptoStim = [];
bl.trialNum = []; bl.usingPanels = [];
for iTrial = 1:numel(bD)
    bl.dffArr(:, :, iTrial) = bD(iTrial).dffMat;                        % --> [volume, glom, trial]
    bl.dffVectAvgRad(:, iTrial) = bD(iTrial).dffVectAvgRad;             % --> [volume, trial]
    bl.dffVectAvgWedge(:, iTrial) = bD(iTrial).dffVectAvgWedge;         % --> [volume, trial]        
    bl.dffVectStrength(:, iTrial) = bD(iTrial).dffVectStrength;         % --> [volume, trial]
    bl.optoStimOnsetTimes{iTrial} = bD(iTrial).optoStimOnsetTimes;      % --> {trial}
    bl.optoStimOffsetTimes{iTrial} = bD(iTrial).optoStimOffsetTimes;   % --> {trial}
    bl.usingOptoStim(iTrial) = bD(iTrial).usingOptoStim;                % --> [trial]
    bl.wedgeDffArr(:, :, iTrial) = bD(iTrial).wedgeDffMat;              % --> [volume, wedge, trial]
    bl.wedgeZscoreArr(:, :, iTrial) = bD(iTrial).wedgeZscoreMat;        % --> [volume, wedge, trial]
    bl.zscoreArr(:, :, iTrial) = bD(iTrial).zscoreMat;                  % --> [volume, glom, trial]
    bl.trialNum(iTrial) = bD(iTrial).trialNum;                          % --> [trial]
    bl.usingPanels(iTrial) = bD(iTrial).usingPanels;                    % --> [trial]
    
    % Trim FicTrac frames if necessary 
    targetFrames = mode(cellfun(@numel, {bD.ftFrameTimes}));
    if numel(bD(iTrial).ftFrameTimes) == targetFrames
       bl.ftFrameTimes = bD(iTrial).ftFrameTimes; 
    end
    % All are --> [frame, trial]
    bl.ftData.intX(:, iTrial) = bD(iTrial).ftData.intX(1:targetFrames);
    bl.ftData.intY(:, iTrial) = bD(iTrial).ftData.intY(1:targetFrames);
    bl.ftData.intHD(:, iTrial) = bD(iTrial).ftData.intHD(1:targetFrames);
    bl.ftData.moveSpeed(:, iTrial) = bD(iTrial).ftData.moveSpeed(1:targetFrames);
    bl.ftData.intFwMove(:, iTrial) = bD(iTrial).ftData.intFwMove(1:targetFrames);
    bl.ftData.intSideMove(:, iTrial) = bD(iTrial).ftData.intSideMove(1:targetFrames);
    bl.ftData.yawSpeed(:, iTrial) = bD(iTrial).ftData.yawSpeed(1:targetFrames);
    bl.ftData.fwSpeed(:, iTrial) = bD(iTrial).ftData.fwSpeed(1:targetFrames);
    bl.ftData.sideSpeed(:, iTrial) = bD(iTrial).ftData.sideSpeed(1:targetFrames);
    
    bl.ftData.moveSpeed(1,iTrial) = bl.ftData.moveSpeed(2,iTrial);
end


% allTrialsWedgeDff = [];
% for iTrial = 1:numel(mD)
%         
%     % Combine with data from other trials
%     allTrialsWedgeDff(:, :, iTrial) = md(iTrial).wedgeDffMat; % --> [volume, wedge, trial]
%     
% end

%% Calculate tuning curves for each wedge across trials

nTrials = numel(bl);

% Get mean panels pos data
panelsPosVols = [];
for iVol = 1:size(bl.wedgeDffArr, 1)
    [~, currVol] = min(abs(bl.panelsFrameTimes - bl.volTimes(iVol)));
    panelsPosVols(iVol) = bl.panelsPosX(currVol);
end
meanDff = [];
for iPos = 1:numel(unique(bl.panelsPosX))
    meanDff(iPos, :, :) = mean(bl.wedgeDffArr(panelsPosVols == (iPos - 1), :, :), 1); % --> [barPos, wedge, trial]    
end
meanDffSmooth = smoothdata(meanDff, 1, 'gaussian', 3);

% Determine subplot grid size
nPlots = size(meanDffSmooth, 2);
if nPlots == 3
    plotPos = [2 2]; % Because [1 3] looks bad
else
    plotPos = numSubplots(nPlots);
end

% Create figure and plots
figure(1);clf;
for iPlot = 1:nPlots
    subplot(plotPos(1), plotPos(2), iPlot);
    
    currData = squeeze(meanDffSmooth(:, iPlot, :)); % --> [barPos, trial]
    currData = currData(:, ~[bl.usingOptoStim]);
    currDataShift = [currData(92:96, :); currData(1:91, :)];
    
    % Plot data
    plotX = -180:3.75:(180 - 3.75);
    imagesc(plotX, [1 size(currData, 2)], ...
            smoothdata(currDataShift, 1, 'gaussian', 3)');
    hold on; plot([0 0], [0 size(currData, 2) + 1], 'color', 'k', 'linewidth', 2)
    ylim([0.5 size(currData, 2) + 0.5])
%     colorbar
    xlabel('Bar position (degrees from front of fly)');
    ylabel('Trial');
    
    % Indicate missing opto stim trials
    skippedTrials = 0;
    optoStimTrials = find([mD.usingOptoStim]);
    for iTrial = 1:numel(optoStimTrials)
        plotY = optoStimTrials(iTrial) - skippedTrials - 0.5;
        skippedTrials = skippedTrials + 1;
        if mD(optoStimTrials(iTrial)).usingPanels
            lineColor = 'green';
        else
            lineColor = 'red';
        end
        xL = xlim;
        plot(xL, [plotY, plotY], 'color', lineColor, 'linewidth', 2);
    end
    
end

%% Plot behavior as a function of bar position


% Get mean panels pos data
panelsPosVols = [];
for iFrame = 1:numel(bl.ftFrameTimes)
    [~, currFrame] = min(abs(bl.panelsFrameTimes - bl.ftFrameTimes(iFrame)));
    panelsPosFrames(currFrame) = bl.panelsPosX(currFrame);
end
meanSpeed = []; meanYawSpeed = []; meanYawVel = [];
for iPos = 1:numel(unique(bl.panelsPosX))
    meanSpeed(iPos, :) = mean(bl.ftData.moveSpeed(panelsPosFrames == (iPos - 1), :), 1);        % --> [barPos, trial]    
    meanYawVel(iPos, :) = mean(bl.ftData.yawSpeed(panelsPosFrames == (iPos - 1), :), 1);        % --> [barPos, trial] 
    meanYawSpeed(iPos, :) = mean(abs(bl.ftData.yawSpeed(panelsPosFrames == (iPos - 1), :)), 1); % --> [barPos, trial]
end


% Shift to make center of plot directly in front of fly
meanSpeed = [meanSpeed(92:96, :); meanSpeed(1:91, :)];
meanYawVel = [meanYawVel(92:96, :); meanYawVel(1:91, :)];
meanYawSpeed = [meanYawSpeed(92:96, :); meanYawSpeed(1:91, :)];

plotVar = meanYawVel;
plotX = -180:3.75:(180 - 3.75);

figure(1);clf;
plotData = plotVar(:, ~[bl.usingOptoStim]);
imagesc(plotX, [1 size(plotData, 2)], ...
    smoothdata(plotData, 1, 'gaussian', 5)');
hold on; plot([0 0], [0 size(plotData, 2) + 1], 'color', 'k', 'linewidth', 2)
ylim([0.5 size(plotData, 2) + 0.5])
colormap('bluewhitered')

% Indicate missing opto stim trials
skippedTrials = 0;
optoStimTrials = find([bl.usingOptoStim]);
for iTrial = 1:numel(optoStimTrials)
    plotY = optoStimTrials(iTrial) - skippedTrials - 0.5;
    skippedTrials = skippedTrials + 1;
    if bl.usingPanels(optoStimTrials(iTrial))
        lineColor = 'green';
    else
        lineColor = 'red';
    end
    xL = xlim;
    plot(xL, [plotY, plotY], 'color', lineColor, 'linewidth', 2);
end



