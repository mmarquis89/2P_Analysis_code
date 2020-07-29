%% Create analysis object
expList = load_expList;
expType = 'PPM2';
expList = expList(contains(expList.expName, expType), :);
a = MoveSpeedAnalysis(expList);

figDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Figs';
if strcmp(expType, 'PPM2')
    neuronType = 'TypeF';
    load(fullfile(figDir, 'allPlotParams_TypeF.mat'), 'allPlotParams');
elseif strcmp(expType, 'D-ANT')
    neuronType = 'TypeD';
    load(fullfile(figDir, 'allPlotParams_TypeD.mat'), 'allPlotParams');
else
    load(fullfile(figDir, 'allPlotParams_TypeD.mat'), 'allPlotParams');
    allPlotParams = allPlotParams(1:size(expList, 1), :);
    allPlotParams.expID = expList.expID;
end
savedAnalysisObjects = [];

%%

expNums = [8];

xLim = [0 300];
yLim = [];
fontSize = 12;
useSavedObj = 0;

% Adjust plot spacing and margins
SV = 0.07;
SH = 0.03;
ML = 0.04;
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
for iExp = 1:numel(expNums)    
    
    % Get current expID and params
    currExpID = allPlotParams.expID{expNums(iExp)};
    disp(currExpID);
    currPlotParams = allPlotParams.params{expNums(iExp)};
    
    currPlotParams.locomotion = 0;
    currPlotParams.odor = 1.5;
    currPlotParams.roiName = 'TypeF-R';
    currPlotParams.ballstop = 0;
    
    a.params = currPlotParams;
    
    % Create axes for current plot
    ax = subaxis(plotLayoutDims(1), plotLayoutDims(2), iExp, ...
            'ml', ML, 'mr', MR, 'mt', MT, 'mb', MB, 'sv', SV, 'sh', SH);
        
    if isempty(savedAnalysisObjects.obj{expNums(iExp)}) || ~useSavedObj
    
        % Initialize filters
        a.filterDefs.expID = currExpID;
        a.filterDefs.roiName = currPlotParams.roiName;
        a = a.init_filters();
        
        % Run analysis
        a = a.analyze();
        
        savedAnalysisObjects.obj{expNums(iExp)} = a;
    else
        % Re-use existing analysis object
        a = savedAnalysisObjects.obj{expNums(iExp)};
    end
  
    a = a.generate_binned_flData('slidingBaseline', 'forward');
    
%     % Plot binned FW speed vs. sliding dF/F
%     a = a.generate_binned_flData('slidingBaseline', 'forward');
%     a.plot_binned_fl(ax);
%     ax.FontSize = fontSize;
%     ax.YLabel.String = 'Mean sliding dF/F';
%     ax.Title.String = [currExpID, ' — ', currPlotParams.roiName];
    
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

    % Yaw vs. abs(forward speed)
    [binMidpoints, binMeans, binSEM] = MoveSpeedAnalysis.bin_data(a.analysisOutput.yawSpeedData, ...
            abs(a.analysisOutput.fwSpeedData), currPlotParams.nAnalysisBins);
    errorbar(ax, binMidpoints, binMeans, binSEM, 'color', 'k', 'linewidth', 1);
    ax.YLabel.String = 'abs(forward speed) (mm/sec)';
    ax.XLabel.String = 'Mean yaw speed (deg/sec)';
    ax.FontSize = fontSize;
    ax.Title.String = [currExpID, ' — ', currPlotParams.roiName];
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

saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Figs';
fileName = ['allExp_binned_yawSpeedVsFwSpeed_plots_', neuronType, '_noLocOnset_scaledX_part_2'];

save_figure(f, saveDir, fileName)









