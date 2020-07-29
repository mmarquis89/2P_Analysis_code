
odorIntegrationPeriodSeconds = 30;
nOdorStimBins = 3;
nSpeedBins = 10;
cm = jet(nOdorStimBins) .* 0.9;
% cm = winter(nOdorStimBins);

legendPos_1 = 'se';
legendPos_2 = 'se';

% Adjust plot spacing and margins
SV = 0.09;
SH = 0.05;
ML = 0.06;
MR = 0.03;
MT = 0.06;
MB = 0.03;

try 
    
% Get source data
speedData = abs(a.analysisOutput.fwSpeedData);
flData = a.analysisOutput.flData;
odorStimData = a.analysisOutput.odorStimData;
volumeRate = a.sourceDataTable.subset.volumeRate(1);
odorIntegrationPeriodVols = round(odorIntegrationPeriodSeconds * volumeRate);

% Calculate odor stim history for each volume
odorStimCounts = [];
odorStimDists = [];
for iVol = 1:numel(flData)
    if iVol > odorIntegrationPeriodVols
        odorStimCounts(iVol) = sum(odorStimData((iVol - odorIntegrationPeriodVols):iVol));
    else
        odorStimCounts(iVol) = sum(odorStimData(1:iVol));
    end
    currOdorStimDist = find(odorStimData(1:iVol), 1, 'last');
    if ~isempty(currOdorStimDist)
        odorStimDists(iVol) = iVol - currOdorStimDist;
    else
        odorStimDists(iVol)= nan;
    end
end

% Create figure
f = figure(19); clf;
f.Color = [1 1 1];
f.UserData.analysisParams.odorIntegrationPeriodSeconds = 100;
f.UserData.analysisParams.nOdorStimBins = 4;
f.UserData.analysisParams.nSpeedBins = 15;
ax = subaxis(4, 2, [1 3], 'ml', ML, 'mr', MR, 'mt', MT, 'mb', MB, 'sv', SV, 'sh', SH);

% Generate histogram bin edges
h = histogram(ax, odorStimCounts, nOdorStimBins);
stimCountBinEdges = h.BinEdges;
h = histogram(ax, odorStimDists, nOdorStimBins);
stimDistBinEdges = h.BinEdges;
h = histogram(ax, speedData, nSpeedBins);
speedBinEdges = h.BinEdges;
cla; hold on;

% Odor count plot
legendStr = {};
for iStimBin = 1:nOdorStimBins
    
    currBinInds = odorStimCounts > stimCountBinEdges(iStimBin) & ...
            odorStimCounts <= stimCountBinEdges(iStimBin + 1);
    currFl = flData(currBinInds);
    currSpeed = speedData(currBinInds);
    legendStr{iStimBin} = [num2str(round(stimCountBinEdges(iStimBin) / volumeRate)), '-', ...
            num2str(round(stimCountBinEdges(iStimBin + 1) / volumeRate)), ...
            ' sec of odor (out of last ', num2str(round(odorIntegrationPeriodVols / ...
            volumeRate)), ')'];
    binMeans = [];
    binMidpoints = [];
    SEM = [];
    for iBin = 1:(numel(speedBinEdges) - 1)
        binMidpoints(iBin) = mean([speedBinEdges(iBin), speedBinEdges(iBin + 1)]);
        currBinData = currFl(currSpeed > speedBinEdges(iBin) & ...
                currSpeed <= speedBinEdges(iBin + 1));
        if numel(currBinData) > 5
            SEM(iBin) = std_err(currBinData);
        else
            SEM(iBin) = nan;
        end
        binMeans(iBin) = mean(currBinData, 'omitnan');
    end
    binMeans(isnan(SEM)) = nan;
    binMidpoints(isnan(SEM)) = nan;
    errorbar(binMidpoints, binMeans, SEM, '-o', 'color', cm(iStimBin, :), 'linewidth', 1)        
end
legend(legendStr, 'location', legendPos_1, 'FontSize', 11)
ax.FontSize = 13;
ylabel('Mean dF/F');
ax.Title.String = {a.filterDefs.expID, 'Grouped by cumulative history of odor exposure'};
xlabel('Speed (mm/sec)');


% Odor distance plot
ax = subaxis(4, 2, [2 4]); hold on;
legendStr = {};
for iStimBin = 1:nOdorStimBins
    
    currBinInds = odorStimDists > stimDistBinEdges(iStimBin) & ...
            odorStimDists <= stimDistBinEdges(iStimBin + 1);
    currFl = flData(currBinInds);
    currSpeed = speedData(currBinInds);
    legendStr{iStimBin} = [num2str(round(stimDistBinEdges(iStimBin) / volumeRate)), '-', ...
            num2str(round(stimDistBinEdges(iStimBin + 1) / volumeRate)), ...
            ' sec since last odor'];
    binMeans = [];
    binMidpoints = [];
    SEM = [];
    for iBin = 1:(numel(speedBinEdges) - 1)
        binMidpoints(iBin) = mean([speedBinEdges(iBin), speedBinEdges(iBin + 1)]);
        currBinData = currFl(currSpeed > speedBinEdges(iBin) & ...
                currSpeed <= speedBinEdges(iBin + 1));
        if numel(currBinData) > 5
            SEM(iBin) = std_err(currBinData);
        else
            SEM(iBin) = nan;
        end
        binMeans(iBin) = mean(currBinData, 'omitnan');
    end
    binMeans(isnan(SEM)) = nan;
    binMidpoints(isnan(SEM)) = nan;
    errorbar(binMidpoints, binMeans, SEM, '-o', 'color', cm(iStimBin, :), 'linewidth', 1)        
end
legend(legendStr, 'location', legendPos_2, 'FontSize', 11)
ax.FontSize = 13;
ax.Title.String = 'Grouped by time since last odor stim';
xlabel('Speed (mm/sec)');

% Visualize integrated odor count values
ax = subaxis(4, 2, 5);
plot(odorStimCounts / volumeRate, 'color', 'b', 'linewidth', 1);
ax.XTickLabel = '';
ax.FontSize = 13;
ylabel('Seconds')
xlabel('Time');

ax = subaxis(4, 2, 7); hold on
plotData = sort(odorStimCounts) / volumeRate; 
plot(plotData, 'color', 'b', 'linewidth', 1.5);
ax.XTickLabel = '';
ax.FontSize = 13;
ylabel('Seconds')
yL = ylim();
for iBin = 1:(nOdorStimBins - 1)
    binLoc = argmin((abs(stimCountBinEdges(iBin + 1) / volumeRate - plotData)));
    plot([binLoc, binLoc], yL, 'color', 'r')
end
title('Red lines = group bin edges', 'FontWeight', 'normal')

% Visualize recent odor distance values
ax = subaxis(4, 2, 6);
plot(odorStimDists / volumeRate, 'color', 'b', 'linewidth', 1);
ax.XTickLabel = '';
ax.FontSize = 13;
xlabel('Time');

ax = subaxis(4, 2, 8); hold on;
plotData = sort(odorStimDists) / volumeRate; 
plot(plotData, 'color', 'b', 'linewidth', 1.5)
ax.XTickLabel = '';
ax.FontSize = 13;
yL = ylim();
for iBin = 1:(nOdorStimBins - 1)
    binLoc = argmin((abs(stimDistBinEdges(iBin + 1) / volumeRate - plotData)));
    plot([binLoc, binLoc], yL, 'color', 'r')
end
title('Red lines = group bin edges', 'FontWeight', 'normal')

catch ME; rethrow(ME); end
    

%% Save figure

saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Figs';
fileName = [a.filterDefs.expID, '_odor_history_dependent_movement_analysis'];

save_figure(f, saveDir, fileName)
