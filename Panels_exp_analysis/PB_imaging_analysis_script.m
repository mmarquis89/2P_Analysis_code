%% LOAD DATA
% Load all analysis data for one or more experiments, and generate arrays of fluorescence data for 
% ordered EB wedges and individual PB glomeruli. 

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData_60D05_7f';
expList = {'20201117-3', '20201117-4', '20201120-2'};
figDir = fullfile(parentDir, 'Figs');

[expMd, trialMd, roiData, ftData, flailingEvents, locEvents, panelsMetadata, wedgeData, glomData] = ...
        load_PB_data(parentDir, expList);
    
% Load file with info about the details of the bath applications
drugTimingMd = readtable(fullfile(parentDir, 'drug_timing_data.csv'), 'delimiter', ',');

%---------------------------------------------------------------------------------------------------

%% CHECK CORRELATION BETWEEN GLOMERULI
% Make sure the glomeruli were identified correctly by looking at the correlation between the
% left and right matching glomeruli.
%---------------------------------------------------------------------------------------------------
expList = unique(expMd.expID);%{'20201210-1'}%

for iExp = 1:numel(expList)
    check_glom_correlation(expList{iExp}, expMd, roiData);
end 

%% PREPARE DATA FOR MULTI-TRIAL ANALYSIS

sourceData = wedgeData; % glomData or wedgeData
flSmWin = 5;

tbl = inner_join(trialMd, expMd, sourceData, panelsMetadata);
tbl = outerjoin(tbl, ftData, 'type', 'left', 'mergekeys', 1);
tbl = outerjoin(tbl, drugTimingMd, 'type', 'left', 'mergekeys', 1);

% Create table of data for bar cycle analysis 
cycleData = [];
for iTrial = 1:size(tbl, 1)
    disp(iTrial);
    if tbl.usingPanels(iTrial)
        currCycleData = get_bar_cycle_data(tbl(iTrial, :), flSmWin);
        cycleData  = [cycleData; currCycleData];
    end
end

% Add fields containing cycle start vols and volTimes to main data table
cycleStartVols = [];
for iTrial = 1:size(tbl, 1)
    if tbl.usingPanels(iTrial)
        currTrialCycles = innerjoin(tbl(iTrial, :), cycleData);
        cycleStartVols{iTrial} = currTrialCycles.cycStartVol;
    else
        cycleStartVols{iTrial} = [];
    end
end
tbl.cycStartVols = cycleStartVols';

% Add a volTimes field to the main data table
volTimes = [];
for iTrial = 1:size(tbl, 1)
    volTimes{iTrial} = (1:tbl.nVolumes(iTrial)) ./ tbl.volumeRate(iTrial);
end
tbl.volTimes = volTimes'; 

%% PLOT OVERVIEW OF VISUAL TUNING AND MOVEMENT FOR A SINGLE TRIAL

expID = '20210118-1';
trialNum = 1;

p = [];
p.flType = 'expDff';  % rawFl, trialDff, or expDff
p.plotPVA = 1;
p.plotMeanPVA = 1;
p.plotBarCycles = 1;
p.useFlow = 1;
p.flMax = [];
p.smWin = 5;
p.figNum = [];

%---------------------------------------------------------------------------------------------------

plot_single_trial_visual_tuning_summary(tbl(strcmp(tbl.expID, expID) & ...
        tbl.trialNum == trialNum, :), p);


%% PLOT SEVERAL STACKED SINGLE-TRIAL HEATMAPS

currExpID = expID;
trialNums = [];

p = [];
p.smWin = 5;
p.flType = 'expDff';
p.plotBarCycles = 1;
p.flMax = [];
p.figPos = [];
p.figNum = 2;
p.SV = 0.002;
p.SH = 0.02;
p.ML = 0.04;
p.MR = 0.02;
p.MT = 0.05;
p.MB = 0.065;

% Separate data for current experiment and trial(s)
if isempty(trialNums)
    trialNums = tbl.trialNum(strcmp(tbl.expID, currExpID));
end
currTbl = tbl(strcmp(tbl.expID, currExpID) & ismember(tbl.trialNum, trialNums), :);

[f, ax] = plot_single_trial_bump_heatmaps(currTbl, p);


%% PLOT TUNING HEATMAPS FOR EACH WEDGE ACROSS TRIALS

currExpID = expID;
trialNums = [];

startTimes = [];
endTimes = [];

p = [];
p.smWin = 5;
p.flType = 'expDff';
p.flMax = [];
p.figPos = [];
p.figNum = 3;
p.SV = 0.01;
p.SH = 0.03;
p.ML = 0.04;
p.MR = 0.02;
p.MT = 0;
p.MB = 0.08;

% Separate data for current experiment and trial(s)
if isempty(trialNums)
    trialNums = tbl.trialNum(strcmp(tbl.expID, currExpID));
end
currTbl = tbl(strcmp(tbl.expID, currExpID) & ismember(tbl.trialNum, trialNums), :);

% Set tuning window start and end times if they were not provided
if isempty(startTimes)
    startTimes = zeros(size(currTbl, 1), 1);
end
if isempty(endTimes)
    endTimes = currTbl.trialDuration; 
end

% Get visual tuning data
timeWindows = table(currTbl.expID, currTbl.trialNum, startTimes, endTimes, 'VariableNames', ...
        {'expID', 'trialNum', 'startTime', 'endTime'});    
tuningTbl = get_mean_visual_tuning(currTbl, timeWindows);

% Condense fluorescence into a single array (dims = [barPos, trial, wedge])
plotFl = reshape(cell2mat(tuningTbl.([p.flType, 'Tuning'])), size(tuningTbl.rawFlTuning{1}, 1), ...
        size(tuningTbl, 1), size(tuningTbl.rawFlTuning{1}, 2));

% Plot data
f = plot_visual_tuning_heatmaps(currTbl, plotFl, p);

% Add title to figure
titleStr = [expID, '  —  EB wedge mean ', p.flType, ' visual tuning'];
if ~isempty(drugTrials)
    titleStr = [titleStr, '  (green lines = ', currTbl.drugName{drugTrials(1)}, ' application)'];
end
h = suptitle(titleStr);
h.FontSize = 18;

%% PLOT ALL TUNING CURVES ON POLAR PLOTS

currExpID = expID;
trialNums = [];

startTimes = [];
endTimes = [];

p = [];
p.smWin = 5;
p.flType = 'expDff';
p.matchRLims = 0;
p.figPos = [];
p.figNum = 4;
p.SV = 0.03;
p.SH = 0.05;
p.ML = 0.03;
p.MR = 0.03;
p.MT = 0.3;
p.MB = 0.05;

% Separate data for current experiment and trial(s)
if isempty(trialNums)
    trialNums = tbl.trialNum(strcmp(tbl.expID, currExpID));
end
currTbl = tbl(strcmp(tbl.expID, currExpID) & ismember(tbl.trialNum, trialNums), :);

% Set tuning window start and end times if they were not provided
if isempty(startTimes)
    startTimes = zeros(size(currTbl, 1), 1);
end
if isempty(endTimes)
    endTimes = currTbl.trialDuration; 
end

% Get visual tuning data
timeWindows = table(currTbl.expID, currTbl.trialNum, startTimes, endTimes, 'VariableNames', ...
        {'expID', 'trialNum', 'startTime', 'endTime'});    
tuningTbl = get_mean_visual_tuning(currTbl, timeWindows);

% Condense fluorescence into a single array (dims = [barPos, trial, wedge])
plotFl = reshape(cell2mat(tuningTbl.([p.flType, 'Tuning'])), size(tuningTbl.rawFlTuning{1}, 1), ...
        size(tuningTbl, 1), size(tuningTbl.rawFlTuning{1}, 2));
    
% Make plots
f = plot_visual_tuning_curves_polar(currTbl, plotFl, p);

% Add title
titleStr = [expID, '  —  EB wedge ', p.flType, ' tuning curves'];
if p.matchRLims
    titleStr = [titleStr, '  (axis limits matched across plots)'];
else
    titleStr = [titleStr, '  (axis limits set independently across plots)'];
end
h = suptitle(titleStr);
h.FontSize = 12;

%% PLOT HEATMAPS ALIGNED BY BAR CYCLE

currExpID = expList{1};
trialNums = [];

p = [];
p.smWin = 5;
p.flType = 'trialDff';
p.flMax = [3];
p.figPos = [];
p.figNum = 5;
p.SV = 0.03;
p.SH = 0.02;
p.ML = 0.03;
p.MR = 0.03;
p.MT = 0.03;
p.MB = 0.05;

% Separate data for current experiment and trial(s)
if isempty(trialNums)
    trialNums = tbl.trialNum(strcmp(tbl.expID, currExpID));
end
currTbl = tbl(strcmp(tbl.expID, currExpID) & ismember(tbl.trialNum, trialNums), :);
currCycleData = cycleData(strcmp(cycleData.expID, currExpID) & ismember(cycleData.trialNum, ...
        trialNums), :); 

f = plot_bar_cycle_heatmaps(currTbl, currCycleData, p);

%% 