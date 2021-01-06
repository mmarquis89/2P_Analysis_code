%% Create analysis object
expList = load_expList;
expType = 'PPM2';
expList = expList(contains(expList.expName, expType), :);

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\new_PPL201_experiments';
expList = {'20201222-1', '20201222-2', '20201228-1', '20201228-2', '20201228-3', '20210102-1', ...
        '20210102-2', '20210102-3', '20210102-4'};

a = MoveSpeedAnalysis(expList', parentDir);

figDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Figs';
if strcmp(expType, 'PPM2')
    neuronType = 'PPM1201';
    load(fullfile(figDir, 'allPlotParams_TypeF.mat'), 'allPlotParams');
elseif strcmp(expType, 'D-ANT')
    neuronType = 'PPM203';
    load(fullfile(figDir, 'allPlotParams_TypeD.mat'), 'allPlotParams');
else
    load(fullfile(figDir, 'allPlotParams_TypeD.mat'), 'allPlotParams');
    allPlotParams = allPlotParams(1:size(expList, 1), :);
    allPlotParams.expID = expList.expID;
end
savedAnalysisObjects = [];

%%

expNums = [];

useSavedObj = 0;
useSavedPlotParams = 0;

speedType = 'move';
roiName = 'TypeD';
flType = 'slidingbaseline'; % 'slidingbaseline', 'normalized', 'trialDff', 'expDff', 'rawFl' 
figDims = [1520 840];

xLim = [0 25];
yLim = [0 0.37];
fontSize = 12;

neuronType = 'PPL203';

% Adjust plot spacing and margins
SV = 0.08;
SH = 0.03;
ML = 0.05;
MR = 0.02;
MT = 0.04;
MB = 0.07;

if isempty(expNums)
   expNums = 1:size(allPlotParams, 1); 
end
if numel(expNums) == 3
    plotLayoutDims = [2, 2];
else
    plotLayoutDims = numSubplots(numel(expNums));
end
if isempty(savedAnalysisObjects) || ~useSavedObj
    if ~useSavedPlotParams
        allPlotParams = table(expList', repmat({[]}, numel(expList), 1), 'variableNames', ...
                {'expID', 'params'});
    end
    savedAnalysisObjects = allPlotParams;
    savedAnalysisObjects.obj = repmat({[]}, size(allPlotParams, 1), 1);
end
f = figure(10); clf;
f.Color = [1 1 1];
if ~isempty(figDims)
    f.Position(3:4) = figDims;
end
allRs = [];
for iExp = 1:numel(expNums)
    
    % Get current expID and params
    currExpID = allPlotParams.expID{expNums(iExp)};
    disp(currExpID);
    if useSavedPlotParams
        currPlotParams = allPlotParams.params{expNums(iExp)};
    else
        currPlotParams = a.params;
    end
    
    currPlotParams.skipTrials = [];
%     if iExp > 1 && iExp < 5
%         currPlotParams.skipTrials(end + 1) = 6;
%     end
    
    currPlotParams.max_moveSpeed = 30;
    currPlotParams.smWinVols = 3;
    currPlotParams.smWinFrames = 5;
    currPlotParams.convertSpeedUnits = 0;
    currPlotParams.locomotion = [1];
    currPlotParams.isolatedmovement = [];
    currPlotParams.grooming = [];
    currPlotParams.odor = 3;
    currPlotParams.quiescence = [];
    currPlotParams.roiName = roiName;
    currPlotParams.ballstop = 0;
    currPlotParams.nHistBins = 50;
    currPlotParams.nAnalysisBins = 25;
    currPlotParams.flType = flType;
    a.params = currPlotParams;
    
    % Create axes for current plot
    ax = subaxis(plotLayoutDims(1), plotLayoutDims(2), iExp, ...
            'ml', ML, 'mr', MR, 'mt', MT, 'mb', MB, 'sv', SV, 'sh', SH);
        
    if isempty(savedAnalysisObjects.obj{expNums(iExp)}) || ~useSavedObj
    
        % Initialize filters
        a.filterDefs.expID = currExpID;
        a.filterDefs.roiName = currPlotParams.roiName;
        a = a.init_filters();
        
        % If the filter matches more than one ROI, just use the first one
        uniqueRois = unique(a.sourceDataTable.subset.roiName);
        if numel(uniqueRois) > 1
            a.filterDefs.roiName = uniqueRois{1};
            a = a.init_filters();
        end
        
        % Run analysis
        a = a.analyze();
        
%         a.plot_speedMat();
%         colormap(viridis)
%         a.plot_flMat();
%         colormap(viridis);
        
        savedAnalysisObjects.obj{expNums(iExp)} = a;
    else
        % Re-use existing analysis object
        a = savedAnalysisObjects.obj{expNums(iExp)};
    end
  %     
    % Plot binned speed vs. dF/F
    a = a.generate_binned_flData(currPlotParams.flType, speedType);
    a.plot_binned_fl(ax);
    ax.FontSize = fontSize;
    ax.YLabel.String = 'Mean dF/F';
    if strcmp(speedType, 'forward')
        cf = a.analysisOutput.r(1);
    else
        cf = a.analysisOutput.r(2);
    end
    allRs(end + 1) = cf;
    ax.Title.String = [currExpID, ' — ', neuronType, ' — r = ', num2str(cf, 2)];
    box off
    
    
%     % Plot binned calculated move speed vs. sliding dF/F
%     a = a.generate_binned_flData('slidingBaseline', 'calculatedMove');
%     a.plot_binned_fl(ax);
%     ax.FontSize = fontSize;
%     ax.YLabel.String = 'Mean sliding dF/F';
%     ax.Title.String = [currExpID, ' — ', currPlotParams.roiName];
    
%     % Plot binned yaw speed vs. sliding dF/F
%     a = a.generate_binned_flData('slidingBaseline', 'yaw');
%     a.plot_binned_fl(ax);
%     ax.FontSize = fontSize;
%     ax.YLabel.String = 'Mean sliding dF/F';
%     ax.Title.String = [currExpID, ' — ', currPlotParams.roiName];

%     % Yaw vs. abs(forward speed)
%     [binMidpoints, binMeans, binSEM] = MoveSpeedAnalysis.bin_data(a.analysisOutput.yawSpeedData, ...
%             abs(a.analysisOutput.fwSpeedData), currPlotParams.nAnalysisBins);
%     errorbar(ax, binMidpoints, binMeans, binSEM, 'color', 'k', 'linewidth', 1);
%     ax.YLabel.String = 'abs(forward speed) (mm/sec)';
%     ax.XLabel.String = 'Mean yaw speed (deg/sec)';
%     ax.FontSize = fontSize;
%     ax.Title.String = [currExpID, ' — ', currPlotParams.roiName];
%     
%     % Forward speed vs. Yaw
%     [binMidpoints, binMeans, binSEM] = MoveSpeedAnalysis.bin_data(a.analysisOutput.fwSpeedData, ...
%             a.analysisOutput.yawSpeedData, currPlotParams.nAnalysisBins);
%     errorbar(ax, binMidpoints, binMeans, binSEM, 'color', 'k', 'linewidth', 1);
%     ax.XLabel.String = 'Forward speed (mm/sec)';
%     ax.YLabel.String = 'Mean yaw speed (deg/sec)';
%     ax.FontSize = fontSize;
%     ax.Title.String = [currExpID, ' — ', currPlotParams.roiName];
    
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
end

%% Save current figure

saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\new_PPL201_experiments\Figs'];
fileName = ['Binned_moveSpeed_vs_slidingBaselineDff_plots_', neuronType, '_matchedAxLims'];
f.UserData.plotParams = currPlotParams;
f.UserData.speedType = speedType;
save_figure(f, saveDir, fileName)









