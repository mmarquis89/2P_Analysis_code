saveFig = 0;

expNums = [];

skipTrials = {11:15, 11:16, 8:15, [1:4, 9, 12:13], 7:10, 8:12, 6:8, 7:10, 8:10, 7:9, 7:11, 7:9, ...
        7:10}';

useSavedObj = 0;

speedType = 'yaw';

plotting_minSpeed = 0;
nAnalysisBins = 15;

roiName = 'FB-DAN';
flType = 'expDff'; % 'slidingbaseline', 'normalized', 'trialDff', 'expDff', 'rawFl' 
figDims = [1750 1000]%[1520 840];
    
xLim = [];
yLim = [];
fontSize = 12;

try 
    
% Adjust plot spacing and margins
SV = 0.08;
SH = 0.03;
ML = 0.05;
MR = 0.02;
MT = 0.04;
MB = 0.07;

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

f = figure(1); clf;
f.Color = [1 1 1];
if ~isempty(figDims)
   f.Position(3:4) = figDims; 
end
allRs = [];
for iExp = 1:numel(expNums)   
    
    % Get current expID and params
    currExpID = expList{expNums(iExp)};
    disp(currExpID);
    
    a.params.roiName = roiName;
    a.params.nHistBins = 25;
    a.params.flType = flType;
    a.params.convertSpeedUnits = 0;
    a.params.smWinVols = 5;
    a.params.smWinFrames = 5;
    a.params.nSmoothRepsFrames = 50;
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

    % Get yawSpeed, fwSpeed, and dF/F data
    yawSpeed = a.analysisOutput.yawSpeedData;
    fwSpeed = a.analysisOutput.fwSpeedData;
    flData = a.analysisOutput.flData;
    
    % Calculate bin edges
    yBinSize = max(yawSpeed) / nAnalysisBins;
    xBinSize = max(fwSpeed) / nAnalysisBins;
    xBinEdges = 0:xBinSize:max(fwSpeed);
    yBinEdges = 0:yBinSize:max(yawSpeed);
    
    
    % Calculate bin averages
    binAverages = nan(nAnalysisBins);
    for xBin = 1:nAnalysisBins
        for yBin = 1:nAnalysisBins
            xBinInds = fwSpeed >= xBinEdges(xBin) & fwSpeed < xBinEdges(xBin + 1);
            yBinInds = yawSpeed >= yBinEdges(yBin) & yawSpeed < yBinEdges(yBin + 1);
            meanVal = mean(flData(xBinInds & yBinInds), 'omitnan');
            if ~isnan(meanVal)
                binAverages(yBin, xBin) = meanVal;
            end
        end
    end
    
    % Plot figure
    imagesc(ax, flipud(binAverages));
%     imagesc(ax, imgaussfilt(flipud(binAverages), 0.5))
%     imagesc(ax, inpaint_nans(flipud(binAverages)))
    colormap(viridis)
    colorbar()
    
    % Format axes
    ax.FontSize = fontSize;
    ax.YLabel.String = 'Yaw speed (deg/sec)';
    ax.XLabel.String = 'Forward speed (mm/sec)';
    if mod(nAnalysisBins, 2)
        midBinInd = ceil(nAnalysisBins / 2);
    else
        midBinInd = nAnalysisBins / 2;
    end
    ax.XTick = [1,  midBinInd - 1, nAnalysisBins];
    ax.XTickLabel = [0, round(mean(xBinEdges((midBinInd - 1):midBinInd))), ...
            round(mean(xBinEdges(end-1:end)))];
    ax.YTick = [1,  midBinInd - 1, nAnalysisBins];
    ax.YTickLabel = [round(mean(yBinEdges(end-1:end))), ...
            round(mean(yBinEdges((midBinInd - 1):midBinInd))), 0];    
   
    % Add titles    
    if strcmp(speedType, 'forward')
        cf = a.analysisOutput.r(1);
    else
        cf = a.analysisOutput.r(2);
    end
    allRs(end + 1) = cf;
    ax.Title.String = [currExpID, ' — ', roiName, ' — r = ', num2str(cf, 2)];
    box off

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

% Set NaN color to black
f.Colormap(1,:) = [0 0 0];
    
% Save figure
if saveFig
    fileName = [roiName, '_2D_heatmap_speed_vs_', flType]
    f.UserData.params = a.params;
    save_figure(f, figDir, fileName);
end

catch ME; rethrow(ME); end