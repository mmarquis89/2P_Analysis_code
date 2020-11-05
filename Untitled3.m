figure(2); clf;
hold on;
bd = a_noPanels.analysisOutput.binnedData;
errorbar(bd.speedBinMidpoints, bd.flBinMeans, bd.flBinSEM, '-o', 'color', 'k', 'linewidth', 1);
bd = a_panels.analysisOutput.binnedData;
errorbar(bd.speedBinMidpoints, bd.flBinMeans, bd.flBinSEM, '-o', 'color', 'r', 'linewidth', 1);

%%

saveFig = 1;

expNums = [1, 2, 5:13];

skipTrials = {11:15, 11:16, 8:15, [1:4, 9, 12:13], 7:10, 8:12, 6:8, 7:10, 8:10, 7:9, 7:11, 7:9, ...
        7:10}';
    
    
panelsSkipTrials = {[1:10, 16:19], [1:10, 14:16],[1:6, 11:16], [1:7, 13:17], 1:5, 1:6, 1:7, 1:6, ...
    1:6, 1:6, 1:6}';    

useSavedObj = 0;

speedType = 'yaw';

plotting_minSpeed = 0;
nAnalysisBins = 45;

roiName = 'EB-DAN';
flType = 'expDff'; % 'slidingbaseline', 'normalized', 'trialDff', 'expDff', 'rawFl' 
figDims = [1750 1000]%[1520 840];
    
xLim = [0 400];
yLim = [];
fontSize = 12;

% Adjust plot spacing and margins
SV = 0.08;
SH = 0.03;
ML = 0.05;
MR = 0.02;
MT = 0.04;
MB = 0.07;

try
if isempty(expNums)
   expNums = 1:numel(expList);
else
    skipTrials = skipTrials(expNums);
end
if numel(expNums) == 3
    plotLayoutDims = [2, 2];
else
    plotLayoutDims = numSubplots(numel(expNums));
end
savedAnalysisObjects = [];
% if isempty(savedAnalysisObjects) || ~useSavedObj
%     savedAnalysisObjects = allPlotParams;
%     savedAnalysisObjects.obj = repmat({[]}, size(allPlotParams, 1), 1);
% end

f = figure(1); clf;
f.Color = [1 1 1];
if ~isempty(figDims)
   f.Position(3:4) = figDims; 
end
allRs = [];
noPanelsBd = {};
panelsBd = {};
for iExp = 1:numel(expNums)   
    
    % Get current expID and params
    currExpID = expList{expNums(iExp)};
    disp(currExpID);
    
    a.params.roiName = roiName;
    a.params.nHistBins = 25;
    a.params.flType = flType;
    a.params.convertSpeedUnits = 0;
    a.params.smWinVols = 5;
    a.params.smWinFrames = 7;
    a.params.nSmoothRepsFrames = 20;
    a.params.skipTrials = skipTrials{iExp};
    a.params.plotting_minSpeed = plotting_minSpeed;
    a.params.nAnalysisBins = nAnalysisBins;
    
    % Create axes for current plot
    ax = subaxis(plotLayoutDims(1), plotLayoutDims(2), iExp, ...
            'ml', ML, 'mr', MR, 'mt', MT, 'mb', MB, 'sv', SV, 'sh', SH);
    
    if ~useSavedObj || isempty(savedAnalysisObjects.obj{expNums(iExp)})
        
        % Initialize filters
        a.filterDefs.expID = currExpID;
        a.filterDefs.roiName = roiName;
        a = a.init_filters();
        
        % If the filter matches more than one ROI, just use the first one
        uniqueRois = unique(a.sourceDataTable.subset.roiName);
        if numel(uniqueRois) > 1
            a.filterDefs.roiName = uniqueRois{1};
            a = a.init_filters();
        end
        
        % Run analysis
        a = a.analyze();
        
        
        savedAnalysisObjects.obj{expNums(iExp)} = a;
    else
        % Re-use existing analysis object
        a = savedAnalysisObjects.obj{expNums(iExp)};
    end
%     
    % Plot binned speed vs. dF/F
    a = a.generate_binned_flData(flType, speedType);
    
    noPanelsBd{iExp} = a.analysisOutput.binnedData;
    a.params.skipTrials = panelsSkipTrials{iExp};
    a = a.analyze();
    a = a.generate_binned_flData(flType, speedType);
    panelsBd{iExp} = a.analysisOutput.binnedData;
    
    % With and without swinging bar stim on panels
    bd = noPanelsBd{iExp};
    errorbar(ax, bd.speedBinMidpoints, bd.flBinMeans, bd.flBinSEM, '-o', 'color', 'k', ...
            'linewidth', 1);
    bd = panelsBd{iExp};
    hold on;
    errorbar(ax, bd.speedBinMidpoints, bd.flBinMeans, bd.flBinSEM, '-o', 'color', 'r', ...
            'linewidth', 1);
    legend('Darkness', 'Swinging bar');
    
%     a.plot_binned_fl(ax);
    ax.FontSize = fontSize;
    ax.YLabel.String = 'Mean dF/F';
    if strcmp(speedType, 'forward')
        cf = a.analysisOutput.r(1);
    else
        cf = a.analysisOutput.r(2);
    end
    allRs(end + 1) = cf;
    ax.Title.String = [currExpID, ' — ', roiName, ' — r = ', num2str(cf, 2)];
    box off
    yStr = flType;
    
%     % Forward speed vs. Yaw
%     [binMidpoints, binMeans, binSEM] = MoveSpeedAnalysis.bin_data(a.analysisOutput.fwSpeedData, ...
%             a.analysisOutput.yawSpeedData, a.params.nAnalysisBins, a.params.plotting_minSpeed);
%     errorbar(ax, binMidpoints, binMeans, binSEM, '-o', 'color', 'k', 'linewidth', 1);
%     ax.XLabel.String = 'Forward speed (mm/sec)';
%     ax.YLabel.String = 'Mean yaw speed (deg/sec)';
%     ax.FontSize = fontSize;
%     ax.Title.String = [currExpID, ' — ', roiName];    
%     yStr = 'yawSpeed';
%     box off
% % 
%      % Yaw vs. forward speed
%     [binMidpoints, binMeans, binSEM] = MoveSpeedAnalysis.bin_data(a.analysisOutput.yawSpeedData, ...
%             a.analysisOutput.fwSpeedData, a.params.nAnalysisBins, a.params.plotting_minSpeed);
%     errorbar(ax, binMidpoints, binMeans, binSEM, '-o', 'color', 'k', 'linewidth', 1);
%     ax.XLabel.String = 'Yaw speed (deg/sec)';
%     ax.YLabel.String = 'Mean forward speed (mm/sec)';
%     ax.FontSize = fontSize;
%     ax.Title.String = [currExpID, ' — ', roiName];   
%     yStr = 'fwSpeed';
%     box off

    % Remove axis labels if not in last row (for X) or first column (for Y)
    if iExp <= numel(expNums) - plotLayoutDims(2)
        ax.XLabel.String = [];
    end
    if ~ismember(iExp, 1:plotLayoutDims(2):numel(expNums))
        ax.YLabel.String = [];
    end

    % Manually set axis limits if necessary
    if ~isempty(xLim)
        xlim(ax, xLim);
    end
    if ~isempty(yLim)
        ylim(ax, yLim);
    end
end%iExp

if saveFig
    if isempty(xLim)
        xlStr = 'scaledX';
    else
        xlStr = 'matchedX';
    end
    if isempty(yLim)
        ylStr = 'scaledY';
    else
        ylStr = 'matchedY';
    end
    fileName = [roiName, '_binned_panels_comparison_', speedType, 'Speed_vs_', yStr, '_', xlStr, '_', ylStr]
    save_figure(f, figDir, fileName);
end


catch ME; rethrow(ME); end