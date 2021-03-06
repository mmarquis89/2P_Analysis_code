%% LOAD DATA
% Load all analysis data for one or more experiments, and generate arrays of fluorescence data for 
% ordered EB wedges and individual PB glomeruli. 

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData_60D05_7f';
% parentDir = 'D:\Dropbox (HMS)\Shared_2p_Data_MY\20210122-1_38A11_P2X2_60D05_7f\ProcessedData';
% expList = {'20201201-1', '20201203-1', '20201203-2', '20201210-1', '20201210-2', '20201201-2', ...
%         '20201117-1', '20201117-3', '20201117-4', '20201120-2', '20201201-3', '20210118-1', ...
%         '20210118-2', '20210119-1', '20201120-1', '20201120-3', '20210122-1', '20210122-2'};
expList = {'20201117-1', '20201117-4', '20201201-1', '20201210-2', '20201120-2'};
expList = sort(expList);

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


sourceData = wedgeData; % glomData or wedgeData
flSmWin = 5;

try 
    
tbl = inner_join(trialMd, expMd, sourceData, panelsMetadata);
tbl = outerjoin(tbl, ftData, 'type', 'left', 'mergekeys', 1);
tbl = outerjoin(tbl, drugTimingMd, 'type', 'left', 'mergekeys', 1);

% Add a volTimes field to the main data table
volTimes = [];
for iTrial = 1:size(tbl, 1)
    volTimes{iTrial} = ((1:tbl.nVolumes(iTrial)) ./ tbl.volumeRate(iTrial))';
end
tbl.volTimes = volTimes';

% Create table of data for bar cycle analysis 
cycleData = [];
pvaOffset = {};
panelsPosVols = {};
bumpAmp = {};
for iTrial = 1:size(tbl, 1)
    currTbl = tbl(iTrial, :);
    disp(['Getting bar cycle data for trial #', num2str(iTrial), ' of ', num2str(size(tbl, 1))]);
    if tbl.usingPanels(iTrial)
        
        % Get bar location for each imaging volume
        panelsPosFrames = currTbl.panelsPosX{:};
        volTimes = tbl.volTimes{iTrial};
        panelsFrameTimes = (1:numel(panelsPosFrames)) ./ currTbl.panelsDisplayRate;
        panelsPosVols{iTrial, 1} = [];
        for iVol = 1:numel(volTimes)
            [~, currVol] = min(abs(panelsFrameTimes - volTimes(iVol)));
            panelsPosVols{iTrial}(iVol, 1) = panelsPosFrames(currVol);
        end
        
        % Get PVA-bar offset data
        pva = currTbl.pvaRad{:} + pi; % Add pi so it ranges from 0-2*pi
        panelsBarPhase = (2*pi * (panelsPosVols{iTrial} ./ max(panelsPosVols{iTrial})));
        smPVA = mod(smoothdata(unwrap(pva), 'gaussian', flSmWin), 2*pi);
        
        % Calculate the bump offset for each volume
        pvaOffset{iTrial, 1} = circ_dist(panelsBarPhase, smPVA);
        
        % Calculate bump amplitude for each volume
        expDff = currTbl.expDff{:};
        bumpAmp{iTrial, 1} = max(expDff, [], 2) - min(expDff, [], 2);
        
        % Get bar cycle data
        currTbl.pvaOffset = {pvaOffset{iTrial}};
        currTbl.bumpAmp = {bumpAmp{iTrial}};
        currCycleData = get_bar_cycle_data(currTbl, flSmWin);
        cycleData  = [cycleData; currCycleData];
            
    else
        pvaOffset{iTrial, 1} = [];
        panelsPosVols{iTrial, 1} = [];
        bumpAmp{iTrial, 1} = [];
    end
end
tbl.pvaOffset = pvaOffset;
tbl.panelsPosVols = panelsPosVols;
tbl.bumpAmp = bumpAmp;

% Add fields containing cycle start vols to the main data table
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

catch ME; rethrow(ME); end

%% PLOT OVERVIEW OF VISUAL TUNING AND MOVEMENT FOR A SINGLE TRIAL

expID = expList{5}
trialNum = 3;

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

f = plot_single_trial_visual_tuning_summary(tbl(strcmp(tbl.expID, expID) & ...
        tbl.trialNum == trialNum, :), p);


%% PLOT SEVERAL STACKED SINGLE-TRIAL HEATMAPS

currExpID = expList{1};
trialNums = [1:3];

saveFig = 0;

p = [];
p.smWin = 5;
p.flType = 'expDff';
p.plotBarCycles = 1;
p.flMax = [];
p.figPos = [1800 950];
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

% Save figure
if saveFig
    f.UserData.plotParams = p;
    figTitle = [currExpID, '_single_trial_heatmaps'];
    save_figure(f, figDir, figTitle);
end

%% PLOT TUNING HEATMAPS FOR EACH WEDGE ACROSS TRIALS

% TODO: add spacer row for darkness trials

currExpID = expList{2};
trialNums = [1 2 4:9];

saveFig = 0;

startTimes = [];
endTimes = [];

p = [];
p.drugTrialStartDelay = 200;
p.smWin = 5;
p.flType = 'expDff';
p.flMax = [];
p.figPos = [1600 950];
p.figNum = 3;
p.SV = 0.01;
p.SH = 0.03;
p.ML = 0.04;
p.MR = 0.02;
p.MT = 0;
p.MB = 0.08;

try 
    
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
drugTrials = find(~isnan(currTbl.startTime));
for iTrial = 1:numel(drugTrials)
    startTimes(drugTrials(iTrial)) = startTimes(drugTrials(iTrial)) + p.drugTrialStartDelay;
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
titleStr = [currExpID, '  �  EB wedge mean ', p.flType, ' visual tuning'];
drugTrials = find(~isnan(currTbl.startTime));
if ~isempty(drugTrials)
    titleStr = [titleStr, '  (green lines = ', currTbl.drugName{drugTrials(1)}, ' application)'];
end
h = suptitle(titleStr);
h.FontSize = 18;

% Save figure
if saveFig
    f.UserData.plotParams = p;
    figTitle = [currExpID, '_visual_tuning_heatmaps'];
    save_figure(f, figDir, figTitle);
end
catch ME; rethrow(ME); end

%% PLOT ALL TUNING CURVES ON POLAR PLOTS


% TODO: add empty axes for visual stim trials

currExpID = expList{1};
trialNums = [1 3 4 5];

saveFig = 1;

startTimes = [];
endTimes = [];

p = [];
p.drugTrialStartDelay = 120;
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

try
    
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
drugTrials = find(~isnan(currTbl.startTime));
for iTrial = 1:numel(drugTrials)
    startTimes(drugTrials(iTrial)) = startTimes(drugTrials(iTrial)) + p.drugTrialStartDelay;
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
titleStr = [currExpID, '  �  EB wedge ', p.flType, ' tuning curves'];
if p.matchRLims
    titleStr = [titleStr, '  (axis limits matched across plots)'];
else
    titleStr = [titleStr, '  (axis limits set independently across plots)'];
end
h = suptitle(titleStr);
h.FontSize = 12;


f.Children(4).Title.String = 'baseline';
f.Children(1).Title.String = 'post-ATP';



% Save figure
if saveFig
    f.UserData.plotParams = p;
    figTitle = [currExpID, '_visual_tuning_polarPlots_second_ATP_trial'];
    save_figure(f, figDir, figTitle);
end
catch ME; rethrow(ME); end

%% PLOT HEATMAPS ALIGNED BY BAR CYCLE

currExpID = expList{3};
trialNums = [];

saveFig = 1;

p = [];
p.smWin = 5;
p.flType = 'expDff';
p.flMax = [];
p.figPos = [1600 950];
p.figNum = 5;
p.SV = 0.03;
p.SH = 0.02;
p.ML = 0.03;
p.MR = 0.03;
p.MT = 0.03;
p.MB = 0.05;

try
    
% Separate data for current experiment and trial(s)
if isempty(trialNums)
    trialNums = tbl.trialNum(strcmp(tbl.expID, currExpID));
end
currTbl = tbl(strcmp(tbl.expID, currExpID) & ismember(tbl.trialNum, trialNums), :);
currCycleData = cycleData(strcmp(cycleData.expID, currExpID) & ismember(cycleData.trialNum, ...
        trialNums), :); 

f = plot_bar_cycle_heatmaps(currTbl, currCycleData, p);

% Add title to figure
titleStr = [currExpID, '  �  EB wedge ', p.flType, ' aligned to bar position'];
drugTrials = find(~isnan(currTbl.startTime));
if ~isempty(drugTrials)
    titleStr = [titleStr, '  (green boxes = ', currTbl.drugName{drugTrials(1)}, ' application, ', ...
            'blue lines = trial bounds)'];
end
h = suptitle(titleStr);
h.FontSize = 18;

% Save figure
if saveFig
    f.UserData.plotParams = p;
    figTitle = [currExpID, '_bar_cycle_aligned_heatmaps'];
    save_figure(f, figDir, figTitle);
end
catch ME; rethrow(ME); end

%% PLOT VECTOR STRENGTH AND OTHER METRICS OF EACH CYCLE THROUGHOUT THE EXPERIMENT

currExpID = expList{3};
trialNums = [];

p = [];
p.smWin = 1;
p.roiNums = [];
p.figPos = [];
p.figNum = 5;
p.SV = 0.03;
p.SH = 0.02;
p.ML = 0.03;
p.MR = 0.03;
p.MT = 0.03;
p.MB = 0.05;

try
    
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

% Get bump amplitude
bumpAmp = [];
for iCycle = 1:size(currCycleData, 1)
    if currCycleData.fullCycle(iCycle)
        bumpAmp(end + 1) = mean(currCycleData.cycBumpAmp{iCycle});
    end
end


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
figure(6); clf
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

% Plot max-min expDff for each bar cycle
figure(8); clf
plot(cycTimes, smoothdata(bumpAmp, 'gaussian', p.smWin), '-o')
hold on;
yL = ylim();
plot([trialStartCycTimes; trialStartCycTimes], repmat(yL', 1, numel(trialStartCycTimes)), ...
        'linewidth', 2, 'color', 'k')
plot_stim_shading([drugStartTimes, drugEndTimes])
xlim([0 cycTimes(end)])
ylim(yL);

catch ME; rethrow(ME); end
