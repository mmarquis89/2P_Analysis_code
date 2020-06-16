%% Create analysis object
expList = load_expList;
expList = expList(contains(expList.expName, 'PPM2'), :);
a = MoveSpeedAnalysis(expList);

%% Set params and run analysis

expNum = 2;

% DataTable filter values
roiName = 'TypeF';
expID = expList.expID{expNum};

% Analysis parameters
smWinVols = 4;
smWinFrames = 6;
nSmoothReps = 15;
lagVols = 2;
flType = 'expDff';
max_moveSpeed = 20;
skipTrials = [];
analysisWinDur = 60;

% Plotting parameters
a.params.nHistBins = 50;
a.params.maxBinCount = 80;
a.params.plotting_minSpeed = 0;
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

a.params.nAnalysisBins = round((numel(a.analysisOutput.flData) / ...
        a.sourceDataTable.subset.volumeRate(1)) / analysisWinDur);
    
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
titleStr = ['dF/F calculated in sliding ', num2str(analysisWinDur), '-sec window'];
title(ax, titleStr);
ax = subaxis(2, 3, 5);
a = a.generate_binned_flData('slidingbaseline');
a.plot_binned_fl(ax);
ax.FontSize = 14;
ax.YLabel.String = 'Mean sliding dF/F';

% Sliding-window normalized Fl and moveSpeed
ax = subaxis(2, 3, 3);
a.plot_2D_hist(ax, 'normalized', 1);
ax.FontSize = 14;
ax.XLabel.String = '';
titleStr = ['Fl and move speed normalized in sliding ', num2str(analysisWinDur), '-sec win'];
title(ax, titleStr);
ax = subaxis(2, 3, 6);
a = a.generate_binned_flData('normalized', 1);
a.plot_binned_fl(ax);
ax.FontSize = 14;

catch ME; rethrow(ME); end

%% Save current figure

saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Figs';
fileName = '';

save_figure(f, saveDir, fileName)




