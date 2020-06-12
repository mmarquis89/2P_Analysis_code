parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';

expList = load_expList();

% Load metadata, FicTrac data, ROI data for all experiments
allExpData = [];
for iExp = 1:size(expList, 1)
    
    currExpID = expList.expID{iExp};
    disp(currExpID);
    
    roiFileName = [currExpID, '_roiData.mat'];
    roiFile = fullfile(parentDir, roiFileName);
    
    ftFileName = [currExpID, '_ficTracData.mat'];
    ftFile = fullfile(parentDir, ftFileName);

    
    expMdFile = fullfile(parentDir, [currExpID, '_expMetadata.csv']);
    trialMdFile = fullfile(parentDir, [currExpID, '_trialMetadata.mat']);

    if exist(roiFile, 'file') && exist(ftFile, 'file') && exist(expMdFile, 'file')
        load(roiFile, 'roiData');
        load(ftFile, 'ftData');
        expMd = readtable(expMdFile, 'delimiter', ',');
        load(trialMdFile, 'trialMetadata');
        trialMd = trialMetadata;
        if exist(odorEventFile, 'file')
            odorEvents = readtable(odorEventFile, 'delimiter', ',');
            currData = inner_join(expMd, trialMd, roiData, ftData);
            allExpData = [allExpData; outerjoin(currData, odorEvents, 'type', 'left', ...
                    'mergekeys', 1)];
        else
            allExpData = [allExpData; inner_join(expMd, trialMd, roiData, ftData)];
            
        end         
    end    
end


% Add odor stim timing where appropriate
odorData = [];
for iExp = 1:size(expList, 1)
    
    currExpID = expList.expID{iExp};
    disp(currExpID);
    
    odorEventFileName = [currExpID, '_event_data_odor.csv'];
    odorEventFile = fullfile(parentDir, odorEventFileName);
    
    if exist(odorEventFile, 'file')
        odorData = [odorData; readtable(odorEventFile, 'delimiter', ',')];
    end    
end

%%

dt = DataTable(allExpData);
filterDefs = struct();
filterDefs.roiName = 'TypeF';
filterDefs.expID = '20190226-3';%[regexclude('0226-3')];   
skipTrials = [2];

max_moveSpeed = 300;
excludeOdorVols = 1;
odorOffsetExcludeTime = 2;
smWinVols = 1;
smWinFrames = 1;
nReps = 15;
lagVols = 2;

nHistBins = 50;
plotting_minSpeed = 0;


dt = dt.initialize_filters(filterDefs);

if skipTrials == 0
   skipTrials = []; 
end

% plots to visualize overall speed data
try
speedMat = cell2mat(dt.subset.moveSpeed')';
speedMat = repeat_smooth(speedMat, nReps, 'Dim', 2, 'smWin', smWinFrames);
speedMat = speedMat * 25 * 4.5;
speedMat(speedMat > max_moveSpeed) = max_moveSpeed;
figure(7); clf; 
imagesc(speedMat)
colorbar()
figure(8);clf; hold on;
plot(sort(speedMat(:)));
plot([0, numel(speedMat)], [1 1] * plotting_minSpeed, 'color', 'r', 'linewidth', 2)
xlim([0 numel(speedMat)])
ylabel('Sorted moveSpeed (mm/sec)');

% Plot to visualize overall raw fluorescence
figure(9);clf;
flMat = cell2mat(dt.subset.rawFl')';
flMat = smoothdata(flMat, 2, 'gaussian', smWinVols);
imagesc(flMat)
colorbar
catch
end

plotX = [];
plotY = [];
subsetData = dt.subset;
if ~isempty(skipTrials)
   subsetData(skipTrials, :) = []; 
end
        
r = []; lags = [];
for iTrial = 1:size(subsetData, 1)
    disp(iTrial)
    currData = subsetData(iTrial, :);

    % Calculate volume times
    volTimes = calc_volTimes(currData.nVolumes, currData.volumeRate, currData.trialDuration, ...
            currData.originalTrialCount);

    % Get smoothed move and yaw speed data
    moveSpeedData = currData.moveSpeed{:};
    moveSpeedData(moveSpeedData == 0) = nan; % Needed for 20181020-1 for some reason
    smMoveSpeed = repeat_smooth(moveSpeedData, nReps, 'Dim', 1, 'smWin', smWinFrames);
    
    % Convert move speed from rad/frame to mm/sec 
    smMoveSpeed = smMoveSpeed * 25 * 4.5;

    smMoveSpeed(smMoveSpeed > max_moveSpeed) = max_moveSpeed;
    
    % Downsample smoothed data to match volTimes    
    frameRate = 1 / median(diff(currData.frameTimes{:}));
    frameTimes = calc_volTimes(numel(currData.frameTimes{:}), frameRate, currData.trialDuration, ...
            currData.originalTrialCount); 
    speedDataVols = [];
    for iVol = 1:numel(volTimes)
        speedDataVols(iVol) = smMoveSpeed(argmin(abs(currData.frameTimes{:} - volTimes(iVol))));
    end
    speedDataVols = smoothdata(speedDataVols', 1, 'gaussian', smWinVols);

    % Get smoothed fl data
    smFl = smoothdata(currData.rawFl{:}, 1, 'gaussian', smWinVols);

    % Drop any nan volumes
    if ~isempty(currData.pmtShutoffVols{:}) && ~any(isnan(currData.pmtShutoffVols{:}))
        speedDataVols(currData.pmtShutoffVols{:}) = nan;
        smFl(currData.pmtShutoffVols{:}) = nan;
    end
    
    % Exclude odor response volumes
    if excludeOdorVols
        currOdorData = innerjoin(currData(:, {'expID', 'trialNum'}), odorData);
        for iStim = 1:size(currOdorData, 1)
            excludeVols = (volTimes > currOdorData(iStim, :).onsetTime & volTimes < ...
                    (currOdorData(iStim, :).offsetTime + odorOffsetExcludeTime));
            speedDataVols(excludeVols) = nan;
            smFl(excludeVols) = nan;
        end        
    end
    
    % Calculate cross-correlation between Fl data and moveSpeed
    [r(iTrial, :), lags(iTrial, :)] = xcorr(smFl(~isnan(smFl+speedDataVols)), speedDataVols(~isnan(smFl+speedDataVols)), 5); 
%     [r(iTrial, :), lags(iTrial, :)] = xcorr(smFl(~isnan(smFl) & speedDataVols > 1), speedDataVols(~isnan(smFl) & speedDataVols > 1), 5);

    % Apply a lag to the speed data
    speedDataVols = speedDataVols(1:end - (lagVols));
    smFl = smFl((lagVols + 1):end);
    
    % Calculate dF/F
    dff = (smFl - currData.trialBaseline) ./ currData.trialBaseline;
    
    plotX = [plotX; speedDataVols];
    plotY = [plotY; dff];
end

% selectVols = 1:4000;
% plotX = plotX(selectVols);
% plotY = plotY(selectVols);

figure(10);clf; hold on;
% plot(lags', r')
plot(lags(1, :), mean(r, 1, 'omitnan'), 'linewidth', 3, 'color', 'k')
set(gca, 'XAxisLocation', 'top')
%
nanVols = isnan(plotX) | isnan(plotY);
plotX(nanVols) = [];
plotY(nanVols) = [];

% Overlay normalized data
figure(1);clf;hold on
set(gcf, 'Color', [1 1 1]);
plot(plotX ./ max(plotX));
plot(plotY ./ max(plotY));
% xlim(xL);

% Scatterplot
figure(2);clf;
set(gcf, 'Color', [1 1 1]);
plot(plotX, plotY, '.')
xlabel('moveSpeed (mm/sec)'); ylabel('dF/F'); 

% 2D histogram
figure(3);clf;
set(gcf, 'Color', [1 1 1]);
plotY(plotX < plotting_minSpeed) = nan;
plotX(plotX < plotting_minSpeed) = nan;
histogram2(plotX, plotY, nHistBins, 'displaystyle', 'tile')
xlabel('moveSpeed (mm/sec)'); ylabel('dF/F'); 

% Separate data into equal-sized moveSpeed bins
plotDataTable = table(plotX, plotY, 'variableNames', {'moveSpeed', 'fl'});
plotDataSort = sortrows(plotDataTable, 'moveSpeed');
binSize = round(size(plotDataSort, 1) / 50);
binEdges = 1:binSize:size(plotDataSort, 1);
flBinMeans = []; SEM = [];
binEdgeMoveSpeeds = [];
for iBin = 1:numel(binEdges)
    if iBin < numel(binEdges)
        binInds = binEdges(iBin):(binEdges(iBin + 1) - 1);
    else 
        binInds = binEdges(iBin):size(plotDataSort, 1);
    end
    binEdgeMoveSpeeds(iBin) = plotDataSort{binInds(end), 'moveSpeed'};
    flBinMeans(iBin) = mean(plotDataSort{binInds, 'fl'}, 'omitnan');
    SEM(iBin) = std_err(plotDataSort{binInds, 'fl'});
end

figure(4);clf; 
plot(binEdgeMoveSpeeds, 'o')
xlabel('Bin number')
ylabel('moveSpeed (mm/sec)')

figure(5);clf; 
plot(binEdgeMoveSpeeds, flBinMeans, '-o')
errorbar(binEdgeMoveSpeeds, flBinMeans, SEM, '-o');
xlabel('binned moveSpeed (mm/sec)')
ylabel('mean bin trial dF/F')

% figure(6);clf; 
% plot(binEdgeMoveSpeeds(1:end-1), flBinMeans(1:end-1), '-o')


