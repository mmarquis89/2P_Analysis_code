% Plot single dF/F heatmap a time window spanning 1+ trials (VISUAL STIM TRIALS)

saveFig = 1;

currExpID = expList{3}
trialNum = 3;
fileNameSuffix = '_PVA_behindFlyBlanked';

p = [];
p.smWin = 5;
p.flType = 'expDff';

try
currTbl = tbl(strcmp(tbl.expID, currExpID) & tbl.trialNum == trialNum, :);
currCyc = cycleData(strcmp(cycleData.expID, currExpID) & cycleData.trialNum == trialNum ...
        & cycleData.fullCycle == 1, :);

% Get Fl data for heatmap
sourceFl = currTbl.(p.flType){:};
smFl = smoothdata(currTbl.(p.flType){:}, 1, 'gaussian', p.smWin);

% Get drug timing data
drugStart = currTbl.startTime;
drugEnd = drugStart + currTbl.duration;

% Get bar position data
barPos = currTbl.panelsPosVols{:};
barPosAngles = [-161.25:3.75:180, -(180-3.75):3.75:-165];
midPos = find(barPosAngles == 0);
barBehindFlyPanelsPositions = [1:8, 81:96] - 1; % Subtracting 1 for zero-indexed
barBehindFlyVols = ismember(currCyc.cycBarPos{1}, barBehindFlyPanelsPositions);
    
% Get mean cycle data (bump amp, PVA offset, vector strength
bumpAmp = cellfun(@(x) mean(x(~barBehindFlyVols)), currCyc.cycBumpAmp);
offset = cellfun(@(x) circ_mean(x(~barBehindFlyVols)), currCyc.cycPvaOffset);

% Get bump offset SD for each cycle
[~, offsetSD] = cellfun(@(x) circ_std(x(~barBehindFlyVols)), currCyc.cycPvaOffset);

% Get vol and cycle times 
volTimes = currTbl.volTimes{:};
cycTimes = cellfun(@mean, currCyc.trialVolTimes);
relVolTimes = volTimes - drugStart;
relCycTimes = cycTimes - drugStart;

% Create and position figure
f = figure(1); clf;
f.Color = [1 1 1];
f.Position(3:4) = [1500 880];
colormap(magma);
clear ax;
nPlots = 5;

% Plot bar position
ax(1) = subaxis(nPlots, 1, 1, 'mb', 0.07, 'ml', 0.06);
plot(relVolTimes, barPosAngles(barPos + 1), 'color', 'k', 'linewidth', 2)
ylabel('barPos');
ax(1).FontSize = 12;
ax(1).XTickLabel = [];
ax(1).YTick = [-180 0 180];
ylim([-180 180])
colorbar();

% Plot dF/F heatmap
ax(2) = subaxis(nPlots, 1, 2);
imagesc([relVolTimes(1) relVolTimes(end)], [1 8], smFl');
ax(2).FontSize = 12;
ax(2).XTick = [];
ax(2).YTick = [];
colorbar();

% Plot PVA
ax(3) = subaxis(nPlots, 1, 3);
barBehindFly = ismember(barPos, barBehindFlyPanelsPositions);
pva = currTbl.pvaRad{:} + pi;
pva = mod(smoothdata(unwrap(pva), 'gaussian', p.smWin), 2*pi);
pva(barBehindFly) = nan;
plot(relVolTimes, rad2deg(pva - pi), 'linewidth', 2);
ylabel('PVA');
ax(3).FontSize = 12;
ax(3).XTickLabel = [];
ylim([-180 180])
colorbar();

% Plot bump amplitude
% ax(3) = subaxis(nPlots, 1, 3);
% plot(relCycTimes, bumpAmp, '-o', 'color', 'b', 'linewidth', 2);
% ylabel('Bump amp');
% ax(3).FontSize = 12;
% ax(3).XTickLabel = [];
% ylim([1 4])
% colorbar();

% Plot offset SD
ax(4) = subaxis(nPlots, 1, 4);
plot(relCycTimes, rad2deg(offsetSD), '-o', 'color', rgb('orange'), 'linewidth', 2);
ylabel('Offset SD (deg)');
ax(4).FontSize = 12;
ax(4).XTickLabel = [];
ax(4).YTick = [30 90 150];
ylim([30 150]);
colorbar();

% Plot mean offset 
ax(5) = subaxis(nPlots, 1, 5);
plot(relCycTimes, rad2deg(offset), '-o', 'color', 'red', 'linewidth', 2);
ylabel('Mean offset (deg)');
ax(5).FontSize = 12;
xlabel('Time (sec)')
ylim([-150 200]);
% ax(5).YTick = 40:20:100;
colorbar();

% Set X lims 
linkaxes(ax, 'x');
xlim([relVolTimes(1), currCyc.trialVolTimes{end}(end-1) - drugStart]);

% Save figure
if saveFig
    f.UserData.expList = expList;
    f.UserData.trialNums = trialNum;
    f.UserData.smWin = p.smWin;
    figTitle = ['EB-DAN_ATP+bar_example_trial_', currExpID, fileNameSuffix];
    save_figure(f, figDir, figTitle);
end

catch ME; rethrow(ME); end




%% Same, but for darkness trials


saveFig = 1;

currExpID = expList{3}
trialNum = 2;

p = [];
p.smWin = 5;
p.flType = 'expDff';
p.bumpAmpAvgWin = 15;

try
currTbl = tbl(strcmp(tbl.expID, currExpID) & tbl.trialNum == trialNum, :);

% Get Fl data for heatmap
sourceFl = currTbl.(p.flType){:};
smFl = smoothdata(currTbl.(p.flType){:}, 1, 'gaussian', p.smWin);

% Get drug timing data
drugStart = currTbl.startTime;
drugEnd = drugStart + currTbl.duration;
    
% Calculate mean bump amp throughout trial
volTimes = currTbl.volTimes{:};
bumpAmp = smoothdata(max(sourceFl, [], 2) - min(sourceFl, [], 2), 1, 'gaussian', p.smWin);
winStartVols = 1:round(p.bumpAmpAvgWin * currTbl.volumeRate):numel(volTimes);
winEndVols = winStartVols(2:end) - 1;
winStartVols = winStartVols(1:end-1);
bumpAmpWin = [];
bumpAmpWinTimes = [];
for iWin = 1:numel(winStartVols)
    currVols = winStartVols(iWin):winEndVols(iWin);
    bumpAmpWin(iWin) = mean(bumpAmp(currVols));
    bumpAmpWinTimes(iWin) = mean(volTimes(currVols));
end

% Get relative volume and bump amp win times
relVolTimes = volTimes - drugStart;
relBumpAmpWinTimes = bumpAmpWinTimes - drugStart;

% Create and position figure
f = figure(1); clf;
f.Color = [1 1 1];
f.Position(3:4) = [1500 880];
colormap(magma);
clear ax;
nPlots = 5; % Keeping the empty spots so I can match the dimensions exactly

% Plot bar position
ax(1) = subaxis(nPlots, 1, 1, 'mb', 0.07, 'ml', 0.06);

% Plot dF/F heatmap
ax(2) = subaxis(nPlots, 1, 2);
imagesc([relVolTimes(1) relVolTimes(end)], [1 8], smFl');
ax(2).FontSize = 12;
ax(2).XTick = [];
ax(2).YTick = [];
colorbar();

% Plot bump amplitude
ax(3) = subaxis(nPlots, 1, 3);
plot(relBumpAmpWinTimes, bumpAmpWin, '-o', 'color', 'b', 'linewidth', 2);
ylabel('Bump amp');
ax(3).FontSize = 12;
ax(3).XTickLabel = [];
ylim([0 5])
colorbar()

% Plot offset SD
ax(4) = subaxis(nPlots, 1, 4);


% Plot mean offset 
ax(5) = subaxis(nPlots, 1, 5);


% Set X lims 
linkaxes(ax, 'x');
xlim([relVolTimes(1), relVolTimes(end)]);

% Save figure
if saveFig
    f.UserData.expList = expList;
    f.UserData.trialNums = trialNum;
    f.UserData.plotParams = p;
    figTitle = ['EB-DAN_ATP+darkness_example_trial_', currExpID];
    save_figure(f, figDir, figTitle);
end

catch ME; rethrow(ME); end




