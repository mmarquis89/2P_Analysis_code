%% Create analysis object
expList = load_expList;
expList = expList(contains(expList.expName, 'D-ANT'), :);
% expList = expList(contains(expList.groupName, 'gapless'), :);
a = MoveSpeedAnalysis(expList);

allPlotParams = [];

%% Set params and run moveSpeed-only analysis

expNum = 61;

% DataTable filter values
roiName = 'TypeD-L';
expID = expList.expID{expNum};

% Analysis parameters
smWinVols = 3;
smWinFrames = 5;
nSmoothReps = 15;
lagVols = 2;
flType = 'rawFl';
                        max_moveSpeed = 35;
                        max_yawSpeed = 1000;
skipTrials = [];
slidingWinDur = 60;
                        nAnalysisBins = 40;

% Plotting parameters
a.params.nHistBins = 50;
                        a.params.maxHistBinCount = 100;
a.params.plotting_minSpeed = -100;
overlayTimeSec = [40];

% Analysis exclusion events
odor = 0.01;
optostim = 1;
soundstim = 1;
grooming = 0;
isolatedmovement = 0;
locomotion = 0;
quiescence = 0;

% Primary analysis plot spacing and margins
SV = 0.05;
SH = 0.05;
ML = 0.05;
MR = 0.015;
MT = 0.06;
MB = 0.07;

try 
disp(expID)

% Set params
a.params.smWinVols = smWinVols;
a.params.smWinFrames = smWinFrames;
a.params.nSmoothReps = nSmoothReps;
a.params.lagVols = lagVols;
a.params.flType = flType;
a.params.max_moveSpeed = max_moveSpeed;
a.params.max_yawSpeed = max_yawSpeed;
a.params.odor = odor;
a.params.optostim = optostim;
a.params.soundstim = soundstim;
a.params.grooming = grooming;
a.params.isolatedmovement = isolatedmovement;
a.params.locomotion = locomotion;
a.params.quiescence = quiescence;
if skipTrials == 0
   skipTrials = []; 
end
a.params.skipTrials = skipTrials;
a.params.nAnalysisBins = nAnalysisBins;
    
% Initialize filters    
a.filterDefs.expID = expID;
a.filterDefs.roiName = roiName;
a = a.init_filters();

% Run analysis
a = a.analyze();

% Plot intial summaries
a.plot_speedMat();
a.plot_sortedSpeed();
a.plot_flMat();

a.sourceDataTable.subset.roiName

a.params.slidingWinDur = slidingWinDur;
a.params.slidingWinSize = round(slidingWinDur * a.sourceDataTable.subset.volumeRate(1));
a.params.nAnalysisBins = nAnalysisBins;

% ------ Plot data from the analysis ------

% Post-analysis summaries
a.xcorr_plot();
a.norm_data_overlay(overlayTimeSec);

% Main analysis plots
a.plot_2D_hist();
a = a.generate_binned_flData();
a.plot_binned_fl();

% Plot several primary output plots in one figure
f = figure(8); clf;
f.Color = [1 1 1];

% Standard dF/F data
ax = subaxis(2, 3, 1, 'ml', ML, 'mr', MR, 'mt', MT, 'mb', MB, 'sv', SV, 'sh', SH);
a.plot_2D_hist(ax);
ax.FontSize = 14;
ax.XLabel.String = '';
titleStr = {[expID, ' — ', roiName, ' — ', a.analysisOutput.paramSnapshot.flType, ' vs. moveSpeed'], ...
        ['2D hist and binned mean ', a.analysisOutput.paramSnapshot.flType]};
title(ax, titleStr);
ax = subaxis(2, 3, 4);
a = a.generate_binned_flData();
a.plot_binned_fl(ax);
ax.FontSize = 14;

% Sliding-window dF/F
ax = subaxis(2, 3, 2);
a.plot_2D_hist(ax, 'slidingbaseline');
ax.FontSize = 14;
ax.XLabel.String = '';
titleStr = ['dF/F calculated in sliding ', num2str(slidingWinDur), '-sec window'];
title(ax, titleStr);
ax = subaxis(2, 3, 5);
a = a.generate_binned_flData('slidingbaseline');
a.plot_binned_fl(ax);
ax.FontSize = 14;
ax.YLabel.String = 'Mean sliding dF/F';

% Sliding-window normalized Fl and moveSpeed
ax = subaxis(2, 3, 3);
a.plot_2D_hist(ax, 'normalized', [], 1);
ax.FontSize = 14;
ax.XLabel.String = '';
titleStr = ['Fl and move speed normalized in sliding ', num2str(slidingWinDur), '-sec win'];
title(ax, titleStr);
ax = subaxis(2, 3, 6);
a = a.generate_binned_flData('normalized', [], 1);
a.plot_binned_fl(ax);
ax.FontSize = 14;

catch ME; rethrow(ME); end

%% MoveSpeed-fwSpeed-yawSpeed analysis comparison

expNum = 24;

% DataTable filter values
roiName = 'TypeD';
expID = expList.expID{expNum};

% Analysis parameters
smWinVols = 3;
smWinFrames = 5;
nSmoothReps = 15;
lagVols = 2;
                        max_moveSpeed = 35;
                        max_yawSpeed = 1000;
skipTrials = [];
slidingWinDur = 60;
                        nAnalysisBins = 15;

% Plotting parameters
a.params.nHistBins = 50;
                        a.params.maxHistBinCount = 100;
a.params.plotting_minSpeed = -100;
overlayTimeSec = [40];

% Analysis exclusion events
odor = 2;
optostim = 1;
soundstim = 1;
grooming = 0;
isolatedmovement = 0;
locomotion = 0;
quiescence = 0;

% Primary analysis plot spacing and margins
SV = 0.05;
SH = 0.05;
ML = 0.04;
MR = 0.02;
MT = 0.07;
MB = 0.08;

try 
disp(expID)

% Set params
a.params.smWinVols = smWinVols;
a.params.smWinFrames = smWinFrames;
a.params.nSmoothReps = nSmoothReps;
a.params.lagVols = lagVols;
a.params.flType = flType;
a.params.max_moveSpeed = max_moveSpeed;
a.params.max_yawSpeed = max_yawSpeed;
a.params.odor = odor;
a.params.optostim = optostim;
a.params.soundstim = soundstim;
a.params.grooming = grooming;
a.params.isolatedmovement = isolatedmovement;
a.params.locomotion = locomotion;
a.params.quiescence = quiescence;
if skipTrials == 0
   skipTrials = []; 
end
a.params.skipTrials = skipTrials;
a.params.nAnalysisBins = nAnalysisBins;
a.params.roiName = roiName;    

% Initialize filters    
a.filterDefs.expID = expID;
a.filterDefs.roiName = roiName;
a = a.init_filters();

% Run analysis
a = a.analyze();

% Plot intial summaries
a.plot_speedMat();
a.plot_sortedSpeed();
a.plot_flMat();

a.sourceDataTable.subset.roiName

a.params.slidingWinDur = slidingWinDur;
a.params.slidingWinSize = round(slidingWinDur * a.sourceDataTable.subset.volumeRate(1));
a.params.nAnalysisBins = nAnalysisBins;

% ------ Plot data from the analysis ------

% Post-analysis summaries
a.xcorr_plot();
a.norm_data_overlay(overlayTimeSec);

% Main analysis plots
a.plot_2D_hist();
a = a.generate_binned_flData();
a.plot_binned_fl();

% Plot several primary output plots in one figure
f = figure(9); clf;
f.Color = [1 1 1];

% Sliding window dF/F vs calculated move speed
ax = subaxis(2, 4, 1, 'ml', ML, 'mr', MR, 'mt', MT, 'mb', MB, 'sv', SV, 'sh', SH);
a.plot_2D_hist(ax, 'slidingbaseline', 'calculatedMove');
ax.FontSize = 14;
ax.XLabel.String = '';
titleStr = {[expID, ' — ', roiName, ' — ', a.analysisOutput.paramSnapshot.flType, ' vs. moveSpeed'], ...
        ['Calc. move speed vs. sliding dF/F']};
title(ax, titleStr);

ax = subaxis(2, 4, 5);
a = a.generate_binned_flData('slidingbaseline', 'calculatedMove');
a.plot_binned_fl(ax);
ax.FontSize = 14;
ax.YLabel.String = 'Mean sliding dF/F';

% Sliding window dF/F vs forward speed
ax = subaxis(2, 4, 2);
a.plot_2D_hist(ax, 'slidingbaseline', 'forward');
ax.FontSize = 14;
ax.XLabel.String = '';
% titleStr = {[expID, ' — ', roiName, ' — ', a.analysisOutput.paramSnapshot.flType, ' vs. moveSpeed'], ...
%         ['2D hist and binned mean ', a.analysisOutput.paramSnapshot.flType]};
titleStr = ['Forward speed vs. sliding dF/F'];
title(ax, titleStr);

ax = subaxis(2, 4, 6);
a = a.generate_binned_flData('slidingbaseline', 'forward');
a.plot_binned_fl(ax);
ax.FontSize = 14;
ax.YLabel.String = 'Mean sliding dF/F';

% Sliding-window dF/F vs Yaw speed
ax = subaxis(2, 4, 3);
a.plot_2D_hist(ax, 'slidingbaseline', 'yaw');
ax.FontSize = 14;
ax.XLabel.String = '';
% titleStr = ['dF/F calculated in sliding ', num2str(slidingWinDur), '-sec window'];
titleStr = ['Yaw speed vs. sliding dF/F'];
title(ax, titleStr);

ax = subaxis(2, 4, 7);
a = a.generate_binned_flData('slidingbaseline', 'yaw');
% a = a.generate_binned_flData('slidingbaseline', 'calculatedMove');
a.plot_binned_fl(ax);
ax.FontSize = 14;
ax.YLabel.String = 'Mean sliding dF/F';

% Forward speed vs. yaw/move speed
ax = subaxis(2, 4, 4);
h = histogram2(a.analysisOutput.fwSpeedData, a.analysisOutput.yawSpeedData, a.params.nHistBins, ...
        'displaystyle', 'tile');
xEdges = h.XBinEdges;
yEdges = h.YBinEdges;
binCounts = h.BinCounts;
if ~isempty(a.params.maxHistBinCount)
    binCounts(binCounts > a.params.maxHistBinCount) = a.params.maxHistBinCount;
end
histogram2(ax, 'XBinEdges', xEdges, 'YBinEdges', yEdges, 'BinCounts', ...
    binCounts, 'displaystyle', 'tile');
ax.XGrid = 'off';
ax.YGrid = 'off';
ylabel('Yaw speed (deg/sec)');
ax.FontSize = 14;
ax.XLabel.String = '';
titleStr = 'Yaw speed vs. forward speed';
title(ax, titleStr);

ax = subaxis(2, 4, 8); hold on
% [binMidpoints, binMeans, binSEM] = MoveSpeedAnalysis.bin_data(a.analysisOutput.fwSpeedData, ...
%         a.analysisOutput.yawSpeedData, nAnalysisBins);
% errorbar(ax, binMidpoints, binMeans, binSEM, 'color', 'k', 'linewidth', 1);
% [binMidpoints, binMeans, binSEM] = MoveSpeedAnalysis.bin_data(a.analysisOutput.moveSpeedData, ...
%         a.analysisOutput.yawSpeedData, nAnalysisBins);
% errorbar(ax, binMidpoints, binMeans, binSEM,'color', 'b', 'linewidth', 1);
[binMidpoints, binMeans, binSEM] = MoveSpeedAnalysis.bin_data(a.analysisOutput.yawSpeedData, ...
        abs(a.analysisOutput.fwSpeedData), nAnalysisBins);
errorbar(ax, binMidpoints, binMeans, binSEM, 'color', 'k', 'linewidth', 1);
[binMidpoints, binMeans, binSEM] = MoveSpeedAnalysis.bin_data(a.analysisOutput.yawSpeedData, ...
        a.analysisOutput.moveSpeedData, nAnalysisBins);
errorbar(ax, binMidpoints, binMeans, binSEM,'color', 'b', 'linewidth', 1);

legend({'abs(forward speed)', 'Move speed'}, 'location', 'nw')

ylabel('Speed (mm/sec)');
xlabel('Mean yaw speed (deg/sec)');

% xlabel('Speed (mm/sec)');
% ylabel('Mean yaw speed (deg/sec)');

% a = a.generate_binned_flData('normalized', [], 1);
% a.plot_binned_fl(ax);
ax.FontSize = 14;

catch ME; rethrow(ME); end

%% Save current figure

saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Figs';
fileName = ['movement_analysis_', roiName, '_noLocOnset_', expID];

f.UserData.analysisParams = a.params;
save_figure(f, saveDir, fileName)

if isempty(allPlotParams)
    allPlotParams = table({expID}, {a.params}, 'variablenames', {'expID', 'params'}); 
else
    allPlotParams = [allPlotParams; table({expID}, {a.params}, 'variablenames', ...
            {'expID', 'params'})]; 
end
