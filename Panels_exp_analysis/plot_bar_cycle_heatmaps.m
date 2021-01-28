function f = plot_bar_cycle_heatmaps(sourceTbl, cycleTbl, plotParams)


src = sourceTbl;
cyc = cycleTbl;
p = plotParams;

drugTrials = find(~isnan(src.startTime));

% Generate array of fl data
trialNums = unique(src.trialNum);
nTrials = numel(trialNums);
maxCycleLen = max(cellfun(@(x) size(x, 1), cyc.trialVolTimes));
trialStartCycles = 1;
plotFl = [];
drugStartCycles = [];
drugEndCycles = [];
cycleCount = 0;
for iTrial = 1:nTrials
    currSrc = src(src.trialNum == trialNums(iTrial), :);
    if currSrc.usingPanels
        currCyc = cyc(cyc.trialNum == trialNums(iTrial), :);
        nCycles = size(currCyc, 1);
        if iTrial < nTrials
            trialStartCycles(end + 1) = trialStartCycles(end) + nCycles;
        end
        for iCycle = 1:nCycles
            cycleCount = cycleCount + 1;
            switch p.flType
                case 'rawFl'
                    currCycFl = currCyc.cycRawFl{iCycle};
                case 'trialDff'
                    currCycFl = currCyc.cycTrialDff{iCycle};
                case 'expDff'
                    currCycFl = currCyc.cycExpDff{iCycle};
            end%Switch
            
            % Pad data with NaNs if necessary
            if size(currCycFl, 1) < maxCycleLen
                currCycFl = [currCycFl; nan(maxCycleLen - size(currCycFl, 1), size(currCycFl, 2))];
            end
            
            % Add to array
            plotFl = cat(3, plotFl, currCycFl); % --> [vol, ROI, cycle]
            
            % Identify drug offset and onsets
            if ismember(iTrial, drugTrials)
                startTime = currSrc.startTime;
                endTime = currSrc.startTime + currSrc.duration;
                currCycVolTimes = currCyc.trialVolTimes{iCycle};
                if startTime > currCycVolTimes(1) & startTime < currCycVolTimes(end)
                    drugStartCycles(end + 1) = cycleCount;
                end
                if endTime > currCycVolTimes(1) & endTime < currCycVolTimes(end)
                    drugEndCycles(end + 1) = cycleCount;
                end
            end
            
        end%iCycle
    else
        fillerData = nan(size(plotFl, 1), size(plotFl, 2), nCycles);
        plotFl = cat(3, plotFl, fillerData);
        if iTrial < nTrials
            trialStartCycles(end + 1) = trialStartCycles(end) + nCycles;
        end
        cycleCount = cycleCount + nCycles;
    end%if
end%iTrial

% Create figure
if isempty(p.figNum)
   p.figNum = 1;
end
f = figure(p.figNum);clf;
if ~isempty(p.figPos)
    f.Position(3:4) = p.figPos;
    if sum(f.Position([2, 4])) > 1080
        f.Position(2) = 50;
    end
end
f.Color = [1 1 1];

% Determine subplot grid size
nPlots = size(plotFl, 2);
if nPlots == 3
    plotPos = [2 2]; % Because [1 3] looks bad
else
    plotPos = numSubplots(nPlots);
end
plotPos = [1 size(plotFl, 2)];

for iRoi = 1:size(plotFl, 2)
    
    ax = subaxis(plotPos(1), plotPos(2), iRoi, 'mt', p.MT, 'mb', p.MB, 'sv', p.SV, 'mr', p.MR, ...
            'ml', p.ML, 'sh', p.SH);
    currFl = smoothdata(squeeze(plotFl(:, iRoi, :)), 1, 'gaussian', p.smWin); % --> [vol, cycle]
    
    if ~isempty(p.flMax)
        currFl(currFl > p.flMax) = p.flMax;
    end
    currFl(currFl < 0) = 0;
    
    barAngles = (-180 + 3.75):3.75:180;
    barAngleLabels = barAngles([5:96, 1:4])
    xx = [barAngles(1), barAngles(end)];
    imagesc(xx, [1:size(currFl, 2)], currFl')
    hold on;
    colormap(magma);
    
    % Plot trial bounds
    xL = xlim();
    for iTrial = 1:numel(trialStartCycles)
        plot(xL, [1, 1] * trialStartCycles(iTrial), 'color', 'c', 'linewidth', 2)
    end
    xlim(xL);
    
    % Indicate transitions to drug application trials
    if ~isempty(drugStartCycles)
        for iDrug = 1:numel(drugStartCycles)
            xL = xlim();
            plot(xL, [-0.5, -0.5] + drugStartCycles(iDrug), 'color', 'g', 'linewidth', 3)
            plot(xL, [0.5, 0.5] + drugEndCycles(iDrug), 'color', 'g', 'linewidth', 3)
            plot([xL(1), xL(1)] + diff(xL)*0.004, [-0.5, + 0.5] + ...
                    [drugStartCycles(iDrug) drugEndCycles(iDrug)], 'color', 'g', 'linewidth', 3)
            plot([xL(2), xL(2)] - diff(xL)*0.004, [-0.5, +0.5] + ...
                    [drugStartCycles(iDrug) drugEndCycles(iDrug)], 'color', 'g', 'linewidth', 3)
        end
    end
    
    ax.YTickLabel = [];
    if iRoi == 1
        ylabel('Bar cycles', 'fontSize', 14)
    end
    ax.XTickLabel = [];
    
end




















end