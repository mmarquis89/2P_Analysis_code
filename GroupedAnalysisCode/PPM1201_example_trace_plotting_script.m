
%% Load data for a list of experiments
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
dataDir = fullfile(parentDir, 'all_experiments');

expIDList = {'20180523-2', '20181020-1', '20190226-3'};

[expMd, trialMd] = load_metadata(expIDList, dataDir);

ftData = load_ft_data(expIDList, dataDir);

roiData = load_roi_data(expIDList, dataDir);

%% Separate data from a target ROI, using the first one if imaging bilaterally

roiName = 'TypeF';

roiNames = regexp(roiData.roiName, [roiName, '-.'], 'match', 'once');
roiNames(cellfun(@isempty, roiNames)) = []; 

tbl = inner_join(expMd, trialMd, ftData, roiData(strcmp(roiData.roiName, roiNames{1}), :));

%%
p = [];

expNum = 3;

p.smWinVols = 5;
p.smWinFrames = 5;
p.smReps = 20;
FRAME_RATE = 25;

p.xLims = [];

% Extract data for current experiment
currExpID = expIDList{expNum};
currExpTbl = tbl(strcmp(tbl.expID, currExpID), :);

% Get all data for the current experiment into single vectors
allFrameTimes = 0;
allVolTimes = 0;
allExpDff = [];
allMoveSpeed = [];
for iTrial = 1:size(currExpTbl, 1)
    currTbl = currExpTbl(currExpTbl.trialNum == iTrial, :);
    allFrameTimes = [allFrameTimes; currTbl.frameTimes{:} + allFrameTimes(end)];
    allVolTimes = [allVolTimes; calc_volTimes(currTbl.nVolumes, currTbl.volumeRate, ...
            currTbl.trialDuration, currTbl.originalTrialCount)' + (currTbl.trialDuration .* ...
            (iTrial - 1))];
    allExpDff = [allExpDff; (smoothdata((currTbl.rawFl{:} - currTbl.expBaseline) ./ ...
            currTbl.expBaseline, 'gaussian', p.smWinVols))];
    allMoveSpeed = [allMoveSpeed; repeat_smooth(currTbl.moveSpeed{:} * 4.5 * FRAME_RATE, p.smReps, ...
            'smWin', p.smWinFrames)];
end
allFrameTimes = allFrameTimes(2:end);
allVolTimes = allVolTimes(2:end);
trialStartVols = 1:currTbl.nVolumes:numel(allVolTimes);

f = figure(6); clf;
f.Color = [1 1 1];

% Plot dF/F
ax_1 = subaxis(2, 1, 1, 'ml', 0.15);
hold on;
plot(allVolTimes, allExpDff, 'linewidth', 2)
yL = ylim;
for iTrial = 1:numel(trialStartVols)
    plot([1, 1] * allVolTimes(trialStartVols(iTrial)), yL, 'color', 'r', 'linewidth', 1);
end
ylim(yL);
if ~isempty(p.xLims)
    xlim(p.xLims)
    ax_1.XTick = [p.xLims(1):2:p.xLims(2)];
    ax_1.XTickLabel = ax_1.XTick - ax_1.XTick(1);
end
ax_1.FontSize = 14;
ylabel('dF/F');
title([currExpTbl.expID{1}, '  —  PPM1201'])

% Plot move speed
ax_2 = subaxis(2, 1, 2);
hold on;
plot(allFrameTimes, allMoveSpeed, 'linewidth', 2)
yL = ylim;
for iTrial = 1:numel(trialStartVols)
    plot([1, 1] * allVolTimes(trialStartVols(iTrial)), yL, 'color', 'r', 'linewidth', 1);
end
ylim(yL);
if ~isempty(p.xLims)
    xlim(p.xLims)
    ax_2.XTick = [p.xLims(1):2:p.xLims(2)];
    ax_2.XTickLabel = ax_2.XTick - ax_2.XTick(1);
end
ylabel('Move speed (mm/sec)')
linkaxes([ax_1, ax_2], 'x')
xlabel('Time (sec)')
ax_2.FontSize = 14;

%% Save figure

figDir = fullfile(parentDir, 'figs');
figTitle = 'PPM1201_speed_modulation_example_trace_4';
f.UserData.plotParams = p;
save_figure(f, figDir, figTitle);

