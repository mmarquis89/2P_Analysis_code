
% Plot all data for each wedge stacked according to bar cycle

currTrial = 3;
prevTrialCycles = 5;
smWin = 5;

f = figure(100); clf;
f.Color = [1 1 1];

nPlots = size(cycleDffData{currTrial}, 2);
if nPlots == 3
    plotPos = [2 2]; % Because [1 3] looks bad
elseif nPlots >= 1
    plotPos = numSubplots(nPlots);
    for iPlot = 1:nPlots
        ax = subaxis(plotPos(1), plotPos(2), iPlot);
        currWedgeData = squeeze(cycleDffData{currTrial}(:, iPlot, :));
        if currTrial > 1 && prevTrialCycles
            prevTrialData = cycleDffData{currTrial - 1};
            if ~isempty(prevTrialData)
                prevWedgeData = squeeze(prevTrialData(:, iPlot, :));
                currWedgeData = [prevWedgeData(:, (end - prevTrialCycles:end)), currWedgeData];
            end
        end
        
        currWedgeData = smoothdata(currWedgeData, 1, 'gaussian', smWin);
        cycleDur = max(cycleVolTimes{currTrial}(end, :) - cycleVolTimes{currTrial}(1, :), ...
                [], 'omitnan');
        imagesc([0, cycleDur], [1, size(currWedgeData, 2)], currWedgeData');
        colormap(magma);
        hold on
        drugStartCycle = floor(currExpData.startTime(currTrial) / cycleDur) + prevTrialCycles;
        relDrugStartTime = mod(currExpData.startTime(currTrial), cycleDur);
        drugEndCycle = floor((currExpData.startTime(currTrial) + currExpData.duration(currTrial)) ...
                 / cycleDur) + prevTrialCycles;
        relDrugEndTime = mod((currExpData.startTime(currTrial) + currExpData.duration(currTrial)), ...
                cycleDur);
        plot([1, 1] * relDrugStartTime, [0, 1] + drugStartCycle - 0.5, 'linewidth', 2, 'color', 'g');
        plot([1, 1] * relDrugEndTime, [0, 1] + drugEndCycle - 0.5, 'linewidth', 2, 'color', 'r');
    end
end

%%


smWin = 5;

allCycleDff = [];
trialStartCycles = 1;
for iTrial = 1:nTrials
    currCycleDff = cycleDffData{iTrial};
    if ~isempty(currCycleDff)
        allCycleDff = cat(3, allCycleDff, currCycleDff);
        trialStartCycles = [trialStartCycles, trialStartCycles(end) + size(currCycleDff, 3)];
    else
        allCycleDff = cat(3, allCycleDff, nan(size(cycleDffData{iTrial - 1})));
        trialStartCycles = [trialStartCycles, trialStartCycles(end) + ...
                size(cycleDffData{iTrial - 1}, 3)];
    end
end

f = figure(100); clf;
f.Color = [1 1 1];

nPlots = size(allCycleDff, 2);
if nPlots == 3
    plotPos = [2 2]; % Because [1 3] looks bad
elseif nPlots >= 1
    plotPos = numSubplots(nPlots);
end
for iPlot = 1:nPlots
    ax = subaxis(plotPos(1), plotPos(2), iPlot);
    currWedgeData = squeeze(allCycleDff(:, iPlot, :));
    currWedgeData = smoothdata(currWedgeData, 1, 'gaussian', smWin);
    imagesc(currWedgeData');
    colormap(magma);
    hold on
    xL = xlim();
    for iTrial = 1:numel(trialStartCycles)
        plot(xL, [1, 1] * trialStartCycles(iTrial), 'color', 'g', 'linewidth', 2)
    end
    xlim(xL);
end

%% Plot vector strength for each bar cycle

smWin = 3;
plotRois = 8;

f = figure(9); clf
hold on;
f.Color = [1 1 1];

cvs = permute(cycleVectorStrength, [2 3 1]); % --> [wedge, cycle, trial] 
cvs = cvs(:, 1:end-1, :);
cvsRs = smoothdata(reshape(cvs, 8, []), 2, 'gaussian', smWin); % --> [wedge, expCycle]
trialBounds = 1:size(cvs, 2):size(cvsRs, 2);
cycleTimes = [];
for iTrial = 1:nTrials
    if currExpData.usingPanels(iTrial)
        cycleStartVolTimes = cycleVolTimes{iTrial}(1, 1:end-1);
    else
        cycleStartVolTimes = nan(1, size(cvs, 2));
    end
    if iTrial > 1
        cycleStartVolTimes = cycleStartVolTimes + sum(currExpData.trialDuration(1:iTrial-1));
    end
    cycleTimes = [cycleTimes, cycleStartVolTimes];
end
for iRoi = plotRois%1:size(test, 1)
    plot(repmat(cycleTimes, 8, 1), cvsRs(iRoi, :), '*', 'color', cm(iRoi, :), 'linewidth', 1)
end
if numel(plotRois) > 1
    plot(repmat(cycleTimes, 8, 1), mean(cvsRs, 1), 'color', 'k', 'linewidth', 3)
end
ax = gca;
ax.XLim = [0, cycleTimes(end)];
for iTrial = 1:numel(trialBounds)
    yL = ylim();
    if strcmp(currExpData.visualStim{iTrial}, 'yes')
        lineColor = 'g';
    else
        lineColor = rgb('yellow');
    end
    if isnan(cycleTimes(trialBounds(iTrial)))
        cycleTimes(trialBounds(iTrial)) = cycleTimes(trialBounds(iTrial) - 1);
    end
    if ismember(iTrial, drugTrials)
        plot([1, 1] * cycleTimes(trialBounds(iTrial)), yL, 'color', lineColor, 'linewidth', 3)
    else
        plot([1, 1] * cycleTimes(trialBounds(iTrial)), yL, 'k', 'linewidth', 3)
    end
    ylim(yL);
end
ax = gca;
ax.XLim = [0, cycleTimes(end)];


f = figure(10); clf
hold on;
f.Color = [1 1 1];
cvp = permute(cycleVectorPhase, [2 3 1]);
cvp = cvp(:, 1:end-1, :);
cvpRs = reshape(cvp, 8, []);
% cvpRs = circ_smooth(reshape(cvp, 8, []), 'dim', 2, 'method', 'gaussian', 'smWin', smWin);
trialBounds = 1:size(cvp, 2):size(cvpRs, 2);
for iRoi = plotRois%1:size(test, 1)
    plotData = cvpRs(iRoi, :);
    plot(repmat(cycleTimes, 8, 1), plotData, '*', 'color', cm(iRoi, :), 'linewidth', 1)
end
ax = gca;
ax.XLim = [0, cycleTimes(end)];
ax.YLim = pi * [0, 2];
for iTrial = 1:numel(trialBounds)
    yL = ylim();
    if strcmp(currExpData.visualStim{iTrial}, 'yes')
        lineColor = 'g';
    else
        lineColor = rgb('yellow');
    end
    if isnan(cycleTimes(trialBounds(iTrial)))
        cycleTimes(trialBounds(iTrial)) = cycleTimes(trialBounds(iTrial) - 1);
    end
    if ismember(iTrial, drugTrials)
        plot([1, 1] * cycleTimes(trialBounds(iTrial)), yL, 'color', lineColor, 'linewidth', 3)
    else
        plot([1, 1] * cycleTimes(trialBounds(iTrial)), yL, 'k', 'linewidth', 3)
    end
    ylim(yL);
end
