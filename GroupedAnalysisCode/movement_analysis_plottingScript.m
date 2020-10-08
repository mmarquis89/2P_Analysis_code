%% Create analysis object
expList = load_expList;
expType = 'PPM2';
expList = expList(contains(expList.expName, expType), :);
a = MoveSpeedAnalysis(expList);

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

expNums = [1:8];

useSavedObj = 1;

speedType = 'move';
roiName = 'TypeF';
flType = 'expDff'; % 'slidingbaseline', 'normalized', 'trialDff', 'expDff', 'rawFl' 
figDims = [];

xLim = [];
yLim = [];
fontSize = 12;

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
    savedAnalysisObjects = allPlotParams;
    savedAnalysisObjects.obj = repmat({[]}, size(allPlotParams, 1), 1);
end
f = figure(1); clf;
f.Color = [1 1 1];
if ~isempty(figDims)
   f.Position(3:4) = figDims; 
end
allRs = [];
for iExp = 1:numel(expNums)    
    
    % Get current expID and params
    currExpID = allPlotParams.expID{expNums(iExp)};
    disp(currExpID);
    currPlotParams = allPlotParams.params{expNums(iExp)};
    
    currPlotParams.locomotion = [];
    currPlotParams.isolatedmovement = [];
    currPlotParams.grooming = [];
    currPlotParams.odor = [];
    currPlotParams.quiescence = [];
    currPlotParams.roiName = roiName;
    currPlotParams.ballstop = 0;
    currPlotParams.nHistBins = 25;
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

saveDir = ['C:\Users\Wilson Lab\Dropbox (Personal)\Marquis ms (LH-DANs)\', ...
        'Manuscript_figure_drafts\Basic_move_speed_analysis'];
fileName = ['Binned_moveSpeed_vs_expDff_plots_', neuronType, '_matchedAxLims'];
f.UserData.plotParams = currPlotParams;
f.UserData.speedType = speedType;
save_figure(f, saveDir, fileName)









