smWin = 3;

trialTbl = tbl(19, :);

% Calculate bump amplitude for each volume
expDff = trialTbl.expDff{:};
bumpAmp = max(expDff, [], 2) - min(expDff, [], 2);

% Get PVA-bar offset data
pva = trialTbl.pvaRad{:} + pi;
panelsPosFrames = trialTbl.panelsPosX{:};
volTimes = (1:numel(pva)) ./ trialTbl.volumeRate;
panelsFrameTimes = (1:numel(panelsPosFrames)) ./ trialTbl.panelsDisplayRate;
panelsPosVols = [];
for iVol = 1:numel(volTimes)
    [~, currVol] = min(abs(panelsFrameTimes - volTimes(iVol)));
    panelsPosVols(iVol) = panelsPosFrames(currVol);
end
panelsBarPhase = (2*pi * (panelsPosVols ./ max(panelsPosVols)));
cycStartVols = trialTbl.cycStartVols{:};
smPVA = mod(smoothdata(unwrap(pva), 'gaussian', smWin), 2*pi);

offset = (circ_dist(panelsBarPhase, smPVA')); + pi;
absOffset = abs(circ_dist(panelsBarPhase, smPVA'));

% barBehindFlyPanelsPositions = [1:8, 81:96] - 1; % Subtract 1 for zero-indexed panels positions
barBehindFlyPanelsPositions = [1:12, 77:96] - 1; % Subtract 1 for zero-indexed panels positions
barBehindFlyVols = ismember(panelsPosVols, barBehindFlyPanelsPositions);
% 
smPVA(barBehindFlyVols) = nan;
offset(barBehindFlyVols) = nan;
bumpAmp(barBehindFlyVols) = nan;

cycMeans = [];
cycStd = [];
cycTimes = [];
cycAmp = [];
for iCyc = 1:numel(cycStartVols)
    if iCyc < numel(cycStartVols)
        currCycVols = cycStartVols(iCyc):(cycStartVols(iCyc + 1) - 1);
        currOffset = offset(currCycVols);
        currOffset = currOffset(~isnan(currOffset));
        cycMeans(iCyc) = circ_mean(currOffset');
        cycStd(iCyc) = circ_std(currOffset');
        cycAmp(iCyc) = mean(bumpAmp(currCycVols(~isnan(currCycVols))), 'omitnan');
        cycTimes(iCyc) = mean(volTimes([cycStartVols(iCyc), cycStartVols(iCyc + 1)]));
    else
        currOffset = offset(cycStartVols(iCyc):end);
        currOffset = currOffset(~isnan(currOffset));
        cycMeans(iCyc) = circ_mean(currOffset');
        cycStd(iCyc) = circ_std(currOffset');
        cycAmp(iCyc) = mean(bumpAmp(cycStartVols(iCyc):end), 'omitnan');
        cycTimes(iCyc) = mean(volTimes([cycStartVols(iCyc), numel(volTimes)]));
    end
end


f = figure(1);clf; 
f.Color = [1 1 1];

% Bar phase and PVA
subaxis(5, 3, 1:3, 'mt', 0.03); hold on;
plot(volTimes, panelsBarPhase, 'color', 'r', 'linewidth', 1);
plot(volTimes, smPVA, 'color', 'k', 'linewidth', 1);
xlim([0 volTimes(end)])
ylim([0 2*pi])
ylabel('Bar and PVA')

% Bump offset
subaxis(5, 3, 4:6); hold on;
plot(volTimes, offset, 'color', 'k', 'linewidth', 1);
plot([0, volTimes(end)], [0 0], '--', 'color', 'b')
% ylim([0 pi]);
yL = ylim();
for iCyc = 1:numel(cycStartVols)
    plot([1 1] * volTimes(cycStartVols(iCyc)), yL, 'color', 'r'); 
end
xlim([0 volTimes(end)])
ylabel('Bump offset')

% Bump amplitude
subaxis(5, 3, 7:9); hold on;
plot(volTimes, bumpAmp, 'color', 'k', 'linewidth', 1);
yL = ylim();
for iCyc = 1:numel(cycStartVols)
    plot([1 1] * volTimes(cycStartVols(iCyc)), yL, 'color', 'r'); 
end
xlim([0 volTimes(end)])
ylim(yL)
ylabel('Bump amp')

% Mean cycle offset
subaxis(5, 3, [10 13]); hold on;
plot(cycTimes, cycMeans, '-o', 'color', 'k', 'linewidth', 3)
xlim([0 volTimes(end)])
% ylim([0 pi]);
if ~isnan(trialTbl.startTime)
    plot([1 1] * trialTbl.startTime, ylim(), 'color', 'b', 'linewidth', 1);
    plot([1 1] * (trialTbl.startTime + trialTbl.duration), ylim(), 'color', 'b', 'linewidth', 1);
end
title('Mean bump offset')


% Cycle offset std_dev
subaxis(5, 3, [11 14]); hold on;
plot(cycTimes, cycStd, '-o', 'color', 'm', 'linewidth', 3)
xlim([0 volTimes(end)])
% ylim([0 1.2]);
if ~isnan(trialTbl.startTime)
    plot([1 1] * trialTbl.startTime, ylim(), 'color', 'b', 'linewidth', 1);
    plot([1 1] * (trialTbl.startTime + trialTbl.duration), ylim(), 'color', 'b', 'linewidth', 1);
end

% Mean cycle bump amplitude
subaxis(5, 3, [12 15]); hold on;
plot(cycTimes, cycAmp, '-o', 'color', rgb('orange'), 'linewidth', 3)
xlim([0 volTimes(end)])
ylim([0 5]);
if ~isnan(trialTbl.startTime)
    plot([1 1] * trialTbl.startTime, ylim(), 'color', 'b', 'linewidth', 1);
    plot([1 1] * (trialTbl.startTime + trialTbl.duration), ylim(), 'color', 'b', 'linewidth', 1);
end
title('Mean bump amplitude')

%% Plot mean 
allExpPanelsPos = cell2mat(tbl.panelsPosVols');
allExpBumpAmp = cell2mat(tbl.bumpAmp);

positions = unique(panelsPosVols);
allExpPosAmps = [];
for iPos = 1:numel(positions)
    allExpPosAmps(iPos) = mean(allExpBumpAmp(allExpPanelsPos == positions(iPos)));
    
end

allExpList = unique(tbl.expID);
expPosAmps = {};
for iExp = 1:numel(allExpList)
   currTbl = tbl(strcmp(tbl.expID, allExpList{iExp}), :);
   currPanelsPos = cell2mat(currTbl.panelsPosVols');
   currBumpAmp = cell2mat(currTbl.bumpAmp);
   currPositions = unique(currPanelsPos);
   
   sortedAmps = sort(currBumpAmp);
   minAmp = sortedAmps(round(numel(sortedAmps)*0.05));
   maxAmp = sortedAmps(round(numel(sortedAmps)*0.95));
   currPosAmps = [];
   expPosAmps{iExp} = [];
   for iPos = 1:numel(currPositions)
       expPosAmps{iExp}(iPos) = mean(currBumpAmp(currPanelsPos == currPositions(iPos)));
       expPosAmpsNorm{iExp}(iPos) = (expPosAmps{iExp}(iPos) - minAmp) ./ maxAmp;
   end
end



startAngle = (-180 + (3.75 * 5));
posAngles = [startAngle:3.75:180, (-180+3.75):3.75:(-180+(4*3.75))];

f = figure(9);clf; hold on;
f.Color = [1 1 1];
plot(posAngles, allExpPosAmps, 'o', 'markersize', 7, 'color', 'k')
yL = ylim();
plot([0 0], yL, '--', 'color', 'r', 'linewidth', 2)
plot([1 1]*-135, yL, '--', 'color', 'g', 'linewidth', 2)
plot([1 1]*135, yL, '--', 'color', 'g', 'linewidth', 2)
ax = gca();
ax.FontSize = 14;
ylabel('Mean bump amplitude (all experiments)')
xlabel('Bar position')
ax.XTick = [-180 -135 -90 -45 0 45 90 135 180];

f = figure(10);clf; hold on;
f.Color = [1 1 1];
for iExp = 1:numel(expPosAmps)
    plot(posAngles, expPosAmps{iExp}, 'o', 'markersize', 7)
end
yL = ylim();
plot([0 0], yL, '--', 'color', 'r', 'linewidth', 1)
plot([1 1]*-135, yL, '--', 'color', 'g', 'linewidth', 1)
plot([1 1]*135, yL, '--', 'color', 'g', 'linewidth', 1)
ax = gca();
ax.FontSize = 14;
ylabel('Mean bump amplitude')
xlabel('Bar position')

f = figure(12);clf;
f.Color = [1 1 1];
ax = subaxis(1, 1, 1, 'mb', 0.05, 'mt', 0.02, 'ml', 0.15); hold on;
expCount = 0;
for iExp = 1:numel(expPosAmps)
    plot(posAngles, expPosAmpsNorm{iExp} + expCount, 'o', 'markersize', 7)
    expCount = expCount + 0.4;
end
yL = [0, expCount + 0.4];
plot([0 0], yL, '--', 'color', 'r', 'linewidth', 1)
plot([1 1]*-135, yL, '--', 'color', 'g', 'linewidth', 1)
plot([1 1]*135, yL, '--', 'color', 'g', 'linewidth', 1)
ylim(yL);
ax = gca();
ax.FontSize = 14;
ylabel('Mean bump amplitude')
xlabel('Bar position')
ax.YTickLabel = [];
ax.XTick = [-135 0 135];

