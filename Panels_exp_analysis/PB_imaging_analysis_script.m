%% LOAD DATA
% Load all analysis data for one or more experiments, and generate arrays of fluorescence data for 
% ordered EB wedges and individual PB glomeruli. 

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData_60D05_7f';
% parentDir = 'D:\Dropbox (HMS)\Shared_2p_Data_MY\20210122-1_38A11_P2X2_60D05_7f\ProcessedData';
expList = {'20210118-2'};
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
expList = unique(expMd.expID);%

for iExp = 1:numel(expList)
    check_glom_correlation(expList{iExp}, expMd, roiData);
end 

%% PREPARE DATA FOR MULTI-TRIAL ANALYSIS

% temp = glomData;
% for i=1:5
%     temp.rawFl{i} = temp.rawFl{i}(:, 1:8);
%     temp.trialDff{i} = temp.trialDff{i}(:, 1:8);
%     temp.expDff{i} = temp.expDff{i}(:, 1:8);
% end
% temp.pvaRad = wedgeData.pvaRad;
% temp.pvaWedge = wedgeData.pvaWedge;

sourceData = wedgeData; % glomData or wedgeData
flSmWin = 5;

tbl = inner_join(trialMd, expMd, sourceData, panelsMetadata);
tbl = outerjoin(tbl, ftData, 'type', 'left', 'mergekeys', 1);
tbl = outerjoin(tbl, drugTimingMd, 'type', 'left', 'mergekeys', 1);

% Create table of data for bar cycle analysis 
cycleData = [];
for iTrial = 1:size(tbl, 1)
    disp(['Getting bar cycle data for trial #', num2str(iTrial), ' of ', num2str(size(tbl, 1))]);
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

expID = expList{1};
trialNum = 2;

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

currExpID = expList{1};
trialNums = [];

p = [];
p.smWin = 5;
p.flType = 'trialDff';
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

currExpID = expList{1};
trialNums = [];

drugTrialStartDelay = 200;

startTimes = [];
endTimes = [];

p = [];
p.smWin = 5;
p.flType = 'trialDff';
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
drugTrials = find(~isnan(tbl.startTime));
for iTrial = 1:numel(drugTrials)
    startTimes(drugTrials(iTrial)) = startTimes(drugTrials(iTrial)) + drugTrialStartDelay;
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
drugTrials = find(~isnan(currTbl.startTime));
if ~isempty(drugTrials)
    titleStr = [titleStr, '  (green lines = ', currTbl.drugName{drugTrials(1)}, ' application)'];
end
h = suptitle(titleStr);
h.FontSize = 18;

%% PLOT ALL TUNING CURVES ON POLAR PLOTS

currExpID = expList{1};
trialNums = [];

drugTrialStartDelay = 200;

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
drugTrials = find(~isnan(tbl.startTime));
for iTrial = 1:numel(drugTrials)
    startTimes(drugTrials(iTrial)) = startTimes(drugTrials(iTrial)) + drugTrialStartDelay;
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
p.flType = 'expDff';
p.flMax = [];
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

%% PLOT VECTOR STRENGTH OF EACH CYCLE THROUGHOUT THE EXPERIMENT

currExpID = expList{1};
trialNums = [];

p = [];
p.smWin = 3;
p.roiNums = [];
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

drugTrials = find(~isnan(currTbl.startTime));

% Get mean cycle times and trial start cycle times for entire experiment
expVolTimes = [];
totalTrialDur = 0;
cycTimes = [];
trialStartCycTimes = [];
for iTrial = 1:size(currTbl, 1)
    currVolTimes = cell2mat(currCycleData(currCycleData.trialNum == trialNums(iTrial), :).trialVolTimes);
    currCycTimes = cellfun(@mean, currCycleData(currCycleData.trialNum == trialNums(iTrial), :).trialVolTimes');
    expVolTimes = [expVolTimes; currVolTimes + totalTrialDur];
    if isempty(currVolTimes)
        trialStartCycTimes(iTrial) = expVolTimes(end);
    else
        trialStartCycTimes(iTrial) = currVolTimes(1) + totalTrialDur;
    end
    cycTimes = [cycTimes, currCycTimes + totalTrialDur];
    totalTrialDur = totalTrialDur + currTbl(currTbl.trialNum == trialNums(iTrial), :).trialDuration;
end

% Get cycle vector and phase data and drop any incomplete cycles
cycVs = cell2mat(currCycleData.cycVectorStrength);
cycVp = cell2mat(currCycleData.cycVectorPhase);
cycVs = cycVs(logical(currCycleData.fullCycle), :);
cycVp = cycVp(logical(currCycleData.fullCycle), :);
cycTimes = cycTimes(logical(currCycleData.fullCycle));

% Get min and max expDff values from each cycle
cycMins = [];
cycMaxes = [];
for iCyc = 1:size(currCycleData, 1)
    if currCycleData.fullCycle(iCyc)
        cycMins(end + 1, :) = min(currCycleData.cycExpDff{iCyc}, [], 'omitnan');
        cycMaxes(end + 1, :) = max(currCycleData.cycExpDff{iCyc}, [], 'omitnan');
    end
end
cycRanges = cycMaxes - cycMins;

% Identify drug application start and end times
drugStartTimes = [];
drugEndTimes = [];
if ~isempty(drugTrials)
    for iTrial = 1:numel(drugTrials)
        drugStartTimes(iTrial, 1) = trialStartCycTimes(drugTrials(iTrial)) + ...
                currTbl(currTbl.trialNum == drugTrials(iTrial), :).startTime;
        drugEndTimes(iTrial, 1) = drugStartTimes(iTrial) + currTbl(currTbl.trialNum == ...
                drugTrials(iTrial), :).duration;
    end
end

% Plot vector strength data
figure(4);clf; 
if isempty(p.roiNums)
   p.roiNums = 1:size(cycVs, 2); 
end
plot(cycTimes, smoothdata(cycVs(:, p.roiNums), 1, 'gaussian', p.smWin), '-o')
hold on;
plot(cycTimes, smoothdata(mean(cycVs(:, p.roiNums), 2), 'gaussian', p.smWin), 'linewidth', 2, 'color', 'k')
yL = ylim();
plot([trialStartCycTimes; trialStartCycTimes], repmat(yL', 1, numel(trialStartCycTimes)), ...
        'linewidth', 2, 'color', 'k')
ylim(yL);
plot_stim_shading([drugStartTimes, drugEndTimes])
xlim([0 cycTimes(end)])


% Plot vector phase data
figure(5); clf
if isempty(p.roiNums)
   p.roiNums = 1:size(cycVp, 2); 
end
plot(cycTimes, cycVp(:, p.roiNums), '*')
hold on;
% plot(cycTimes, smoothdata(mean(cycVs(:, p.roiNums), 2), 'gaussian', p.smWin), 'linewidth', 2, 'color', 'k')
yL = pi * [0 2];
ylim(yL);
plot([trialStartCycTimes; trialStartCycTimes], repmat(yL', 1, numel(trialStartCycTimes)), ...
        'linewidth', 2, 'color', 'k')
plot_stim_shading([drugStartTimes, drugEndTimes])
xlim([0 cycTimes(end)])

% Plot max-min expDff for each bar cycle
figure(5); clf
if isempty(p.roiNums)
   p.roiNums = 1:size(cycRanges, 2); 
end
plot(cycTimes, cycRanges(:, p.roiNums), '*')
hold on;
% plot(cycTimes, smoothdata(mean(cycVs(:, p.roiNums), 2), 'gaussian', p.smWin), 'linewidth', 2, 'color', 'k')
yL = ylim();
plot([trialStartCycTimes; trialStartCycTimes], repmat(yL', 1, numel(trialStartCycTimes)), ...
        'linewidth', 2, 'color', 'k')
plot_stim_shading([drugStartTimes, drugEndTimes])
xlim([0 cycTimes(end)])
ylim(yL);



