
% Plot all data for each wedge stacked according to bar cycle

currTrial = 8;
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

smWin = 1;
plotRois = 1:8;

f = figure(9); clf
hold on;
f.Color = [1 1 1];

test = permute(cycleVectorStrength, [2 3 1]);
test = test(:, 1:end-1, :);
test2 = smoothdata(reshape(test, 8, []), 2, 'gaussian', smWin);
trialBounds = 1:size(test, 2):size(test2, 2);
for iRoi = plotRois%1:size(test, 1)
    plot(cycleTimes, test2(iRoi, :), '-o', 'color', cm(iRoi, :), 'linewidth', 1)
end
if numel(plotRois) > 1
    plot(cycleTimes, mean(test2, 1), 'color', 'k', 'linewidth', 5)
end
% for iTrial = 1:numel(trialBounds)
%     yL = ylim();
%     if strcmp(currExpData.visualStim{iTrial}, 'yes')
%         lineColor = 'g';
%     else
%         lineColor = rgb('yellow');
%     end
%     if ismember(iTrial, drugTrials)
%         plot([1, 1] * trialBounds(iTrial), yL, 'color', lineColor, 'linewidth', 3)
%     else
%         plot([1, 1] * trialBounds(iTrial), yL, 'k', 'linewidth', 3)
%     end
%     ylim(yL);
% end

f = figure(10); clf
hold on;
f.Color = [1 1 1];

test = permute(cycleVectorPhase, [2 3 1]);
test = test(:, 1:end-1, :);
test2 = smoothdata(reshape(test, 8, []), 2, 'gaussian', smWin);
trialBounds = 1:size(test, 2):size(test2, 2);
for iRoi = plotRois%1:size(test, 1)
    plotData = test2(iRoi, :);
    plotData(plotData > pi) = plotData(plotData > pi) - pi;
    plot((2 * plotData), '-o', 'color', cm(iRoi, :), 'linewidth', 1)
end
% plot(mean(test2, 1), 'color', 'k', 'linewidth', 5)
for iTrial = 1:numel(trialBounds)
    yL = ylim();
    if strcmp(currExpData.visualStim{iTrial}, 'yes')
        lineColor = 'g';
    else
        lineColor = rgb('yellow');
    end
    if ismember(iTrial, drugTrials)
        plot([1, 1] * trialBounds(iTrial), yL, 'color', lineColor, 'linewidth', 3)
    else
        plot([1, 1] * trialBounds(iTrial), yL, 'k', 'linewidth', 3)
    end
    ylim(yL);
end


