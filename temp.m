

% sourceData = bl.wedgeRawFlArr;
sourceData = bl.wedgeDffArr;
% sourceData = bl.wedgeZscoreArr;

showStimTrials = 0;


nTrials = numel(bl);

% Get mean panels pos data
panelsPosVols = [];
for iVol = 1:size(sourceData, 1)
    [~, currVol] = min(abs(bl.panelsFrameTimes - bl.volTimes(iVol)));
    panelsPosVols(iVol) = bl.panelsPosX(currVol);
end
meanData = [];
for iPos = 1:numel(unique(bl.panelsPosX))
    meanData(iPos, :, :) = ...
            mean(sourceData(panelsPosVols == (iPos - 1), :, :), 1); % --> [barPos, wedge, trial]    
end
meanDataSmooth = smoothdata(meanData, 1, 'gaussian', 3);
meanDataShift = cat(1, meanDataSmooth(92:96, :, :), meanDataSmooth(1:91, :, :));

plotX = -180:3.75:(180 - 3.75);

%%

currTrial = 1;



currData = meanDataShift(:, :, currTrial);

f = figure(2);clf;
f.Color = [1 1 1];
imagesc(plotX, [1 size(currData, 2)], ...
        smoothdata(currData, 1, 'gaussian', 3)');
hold on; plot([0 0], [0 size(currData, 2) + 1], 'color', 'k', 'linewidth', 2)
ylim([0.5 size(currData, 2) + 0.5])

















