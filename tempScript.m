%% Create analysis object
expList = load_expList;
expList = expList(contains(expList.expName, 'PPM2'), :);
a = MoveSpeedAnalysis(expList);

%% Set params and run analysis

expNum = 10;

% DataTable filter values
roiName = 'TypeF';
expID = expList.expID{expNum};
disp(expID)

% Analysis parameters
smWinVols = 1;
smWinFrames = 1;
nSmoothReps = 15;
lagVols = 2;
flType = 'expDff';
max_moveSpeed = 30;
skipTrials = [2];

% Analysis exclusion events
odor = 2;
optostim = 1;
soundstim = 1;
grooming = 0;
isolatedmovement = 0;
locomotion = 0;
quiescence = 0;

try 

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

% Post-analysis summaries
a.xcorr_plot();
a.norm_data_overlay(nOverlaySamples);

% Main analysis plots
a.plot_2D_hist();
a = a.generate_binned_flData();
a.plot_binned_fl();

catch ME; rethrow(ME); end

%% Plot data from the analysis

% Plotting parameters
a.params.nHistBins = 50;
a.params.nAnalysisBins = 30;
a.params.plotting_minSpeed = 0;
nOverlaySamples = [1000];

% Plot intial summaries
a.params.flType = 'expDff';
a.plot_speedMat();
a.plot_sortedSpeed();
a.plot_flMat();

% Post-analysis summaries
a.xcorr_plot();
a.norm_data_overlay(nOverlaySamples);

% Main analysis plots
a.plot_2D_hist();
a = a.generate_binned_flData();
a.plot_binned_fl();






