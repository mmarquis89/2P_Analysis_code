%% Extract single-cycle summary data for all experiments

% 5-5-5-5-5 timing
darkExpList = {'20201201-1', '20201203-1', '20201203-2', '20201210-1', '20201210-2', '20201201-2'};
darkExpTrialNums = {3:7, 3:7, 1:5, 3:7, 1:5, 2:6};

% % 5-5-10 timing
% visExpList = {'20201117-1', '20201117-3', '20201120-2', '20210118-1', '20210119-1'};
% visExpTrialNums = {1:4, 2:5, 2:5, 2:5, 2:6};

% 5-5-5 timing
visExpList = {'20201117-1', '20201117-3', '20201117-4', '20201120-2', '20201201-3', '20210118-1', ...
    '20210118-2', '20210119-1'};
visExpTrialNums = {1:3, 2:4, 2:4, 2:4, 3:5, 2:4, 2:4, 2:4};

% % 10-5-5 timing
% visExpList = {'20201117-3', '20201117-4', '20201120-2', '20210118-1', ...
%     '20210118-2', '20210119-1'};
% visExpTrialNums = {1:4, 1:4, 1:4, 1:4, 1:4, 1:4};

% % 10-5-10 timing
% visExpList = {'20201117-3', '20201120-2', '20210118-1', '20210119-1'};
% visExpTrialNums = {1:5, 1:5, 1:5, 1:5};

% % "First vis" 5-5-5 timing
% visExpList = {'20201117-1', '20201117-3', '20201117-4', '20201120-2', '20201201-3', '20210118-1', ...
%     '20210118-2', '20210119-1', '20201201-1', '20201203-1', '20201203-2', '20201210-1', ...
%     '20201210-2', '20201201-2', '20210122-1', '20210122-2'};
% visExpTrialNums = {1:3, 2:4, 2:4, 2:4, 3:5, 2:4, 2:4, 2:4, 5:7, 5:7, 3:5, 5:7, 3:5, 4:6, 6:8, 5:7};

% % "Every possible vis" 5-5-5 timing
% visExpList = {  '20201117-1', '20201117-3', '20201117-4', '20201120-2', '20201201-3', ...
%                 '20210118-1', '20210118-2', '20210119-1', '20201201-1', '20201203-1', ...
%                 '20201203-2', '20201210-1', '20201210-2', '20201201-2', '20210122-1', ...
%                 '20210122-2', '20201117-1', '20201117-3', '20201117-4', '20201120-2', ...
%                 '20210118-1', '20201203-2'};
% visExpTrialNums = {1:3, 2:4, 2:4, 2:4, 3:5, ...
%                    2:4, 2:4, 2:4, 5:7, 5:7, ...
%                    3:5, 5:7, 3:5, 4:6, 6:8, ...
%                    5:7, 4:6, 5:7, 4:6, 5:7, ...
%                    5:7, 7:9};
% 
% % 10-5-5 vis genetic control
% visExpList = {'20201120-1', '20201120-3'};
% visExpTrialNums = {1:4, 1:4};

visExpTrialNums = visExpTrialNums(~strcmp(visExpList, '20210118-2'));
visExpList = visExpList(~strcmp(visExpList, '20210118-2'));
% 
% currExpList = visExpList;
% trialNums = visExpTrialNums;
currExpList = darkExpList;
trialNums = darkExpTrialNums;

try
allCycRoiData = [];
for iExp = 1:numel(currExpList)
    
    currExpID = currExpList{iExp};
    currTrialNums = trialNums{iExp};   

    % Separate data for current experiment and trial(s)
    currTbl = tbl(strcmp(tbl.expID, currExpID) & ismember(tbl.trialNum, currTrialNums), :);
    currCycleData = cycleData(strcmp(cycleData.expID, currExpID) & ismember(cycleData.trialNum, ...
        currTrialNums), :);
    
    % Get cycle vector and phase data and drop any incomplete cycles
    cycVs = cell2mat(currCycleData.cycVectorStrength);
    cycVp = cell2mat(currCycleData.cycVectorPhase);
    cycVs = cycVs(logical(currCycleData.fullCycle), :);
    cycVp = cycVp(logical(currCycleData.fullCycle), :);

    
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

    % Create new row for each ROI and append to the appropriate table
    for iRoi = 1:size(cycVp, 2)
        newRow = table({currExpID}, iRoi, 'variableNames', {'expID', 'roiNum'});
        newRow.cycVs = {cycVs(:, iRoi)};
        newRow.cycVp = {cycVp(:, iRoi)};
        newRow.cycMin = {cycMins(:, iRoi)};
        newRow.cycMax = {cycMaxes(:, iRoi)};
        newRow.cycRange = {cycRanges(:, iRoi)};
        allCycRoiData = [allCycRoiData; newRow];
    end
    allCycRoiData = [allCycRoiData; newRow];
end


allCycExpData = [];
for iExp = 1:numel(currExpList)
    
    currExpID = currExpList{iExp};
    currTrialNums = trialNums{iExp};   

    % Separate data for current experiment and trial(s)
    currTbl = tbl(strcmp(tbl.expID, currExpID) & ismember(tbl.trialNum, currTrialNums), :);
    currCycleData = cycleData(strcmp(cycleData.expID, currExpID) & ismember(cycleData.trialNum, ...
        currTrialNums), :);
    
    % Get mean cycle times and trial start cycle times for entire experiment
    expVolTimes = [];
    totalTrialDur = 0;
    cycTimes = [];
    trialStartCycTimes = [];
    for iTrial = 1:size(currTbl, 1)
        currVolTimes = cell2mat(currCycleData(currCycleData.trialNum == currTrialNums(iTrial), :).trialVolTimes);
        currCycTimes = cellfun(@mean, currCycleData(currCycleData.trialNum == currTrialNums(iTrial), :).trialVolTimes');
        expVolTimes = [expVolTimes; currVolTimes + totalTrialDur];
        if isempty(currVolTimes)
            trialStartCycTimes(iTrial) = expVolTimes(end);
        else
            trialStartCycTimes(iTrial) = currVolTimes(1) + totalTrialDur;
        end
        cycTimes = [cycTimes, currCycTimes + totalTrialDur];
        totalTrialDur = totalTrialDur + currTbl(currTbl.trialNum == currTrialNums(iTrial), :).trialDuration;
    end
    cycTimes = cycTimes(logical(currCycleData.fullCycle));
    
    % Identify drug application start and end times
    drugTrials = find(~isnan(currTbl.startTime));
    drugStartTimes = [];
    drugEndTimes = [];
    if ~isempty(drugTrials)
        for iTrial = 1:numel(drugTrials)
            drugStartTimes(iTrial, 1) = trialStartCycTimes(drugTrials(iTrial)) + ...
                currTbl(drugTrials(iTrial), :).startTime;
            drugEndTimes(iTrial, 1) = drugStartTimes(iTrial) + currTbl(drugTrials(iTrial), ...
                    :).duration;
        end
    end
    
    % Get mean bump amplitude and mean/STD bar phase-PVA offset for each cycle
    cycBumpAmp = [];
    cycOffset = [];
    cycOffsetStd = [];
    for iCyc = 1:size(currCycleData, 1)
        if currCycleData.fullCycle(iCyc)
            
            % Identify bar positions that are behind the fly
            barBehindFlyPanelsPositions = [1:8, 81:96] - 1; % Subtract 1 for zero-indexed panels positions
            %     barBehindFlyPanelsPositions = [1:12, 77:96] - 1; % Subtract 1 for zero-indexed panels positions
            barBehindFlyVols = ismember(currCycleData.cycBarPos{iCyc}, barBehindFlyPanelsPositions);
            
            cycOffset(end + 1, :) = ...
                    circ_mean(currCycleData.cycPvaOffset{iCyc}(~logical(barBehindFlyVols)));
            [~, cycOffsetStd(end + 1, :)] = ...
                    circ_std(currCycleData.cycPvaOffset{iCyc}(~logical(barBehindFlyVols)) + pi);
            cycBumpAmp(end + 1, :) = ...
                    mean(currCycleData.cycBumpAmp{iCyc}(~logical(barBehindFlyVols)));              
        end%if
    end%iCyc
    
    newRow = table({currExpID}, 'variableNames', {'expID'});
    newRow.cycTimes = {cycTimes'};
    newRow.drugStartTimes = {drugStartTimes};
    newRow.drugEndTimes = {drugEndTimes};
    newRow.cycOffset = {cycOffset};
    newRow.cycOffsetStd = {cycOffsetStd};
    newRow.cycBumpAmp = {cycBumpAmp};
    allCycExpData = [allCycExpData; newRow];
end%iExp

catch ME; rethrow(ME); end



%%

smWin = 3;
    singleRois = 0;
singleExpts = 1;
shadeSEM = 1;
shadeExpSEM = 0;

try
cyc = innerjoin(allCycRoiData, allCycExpData);

cycTimes = cyc.cycTimes{1};
Vs = cell2mat(cyc.cycVs');
Vp = cell2mat(cyc.cycVp');
cycMin = cell2mat(cyc.cycMin');
cycMax = cell2mat(cyc.cycMax');
cycRanges = cell2mat(cyc.cycRange');

% Vector strength
f = figure(1); clf; hold on;
f.Color = [1 1 1];
cm = lines(numel(unique(cyc.expID)));
allExpIDs = unique(cyc.expID);
xx = cycTimes - cyc.drugStartTimes{1}(1);
legendStr = {};
for iExp = 1:numel(allExpIDs)
    currExpID = allExpIDs{iExp};
    currCyc = cyc(strcmp(cyc.expID, currExpID), :);
    if singleRois
        for iRoi = 1:size(currCyc, 1)
            plot(xx, smoothdata(currCyc.cycVs{iRoi}, 1, 'gaussian', smWin), 'color', ...
                    cm(iExp, :));
        end
    elseif singleExpts
        currData = smoothdata(cell2mat(currCyc.cycVs'), 1, 'gaussian', smWin);
        SE = std_err(currData, 2);
        plot(xx, mean(currData, 2), 'linewidth', 2, 'color', cm(iExp, :));
        legendStr{end + 1} = currCyc.expID{1};
        if shadeExpSEM
            upperY = (mean(currData, 2) + SE);
            lowerY = (mean(currData, 2) - SE);
            jbfill(xx', upperY', lowerY', cm(iExp, :), cm(iExp, :), 1, 0.2);
        end
    end
end
legend(legendStr, 'autoupdate', 'off');
meanData = mean(smoothdata(Vs, 1, 'gaussian', smWin), 2);
plot(xx, meanData, 'color', 'k', 'linewidth', 3);
if shadeSEM
    SE = std_err(smoothdata(Vs, 1, 'gaussian', smWin), 2);
    upperY = (meanData + SE);
    lowerY = (meanData - SE);
    jbfill(xx', upperY', lowerY', rgb('black'), rgb('black'), 1, 0.2);
end
if numel(cyc.drugStartTimes{1}) > 1
    shadeColor = rgb('orange');
else
    shadeColor = rgb('green');
end
plot_stim_shading([0, cyc.drugEndTimes{1}(1) - cyc.drugStartTimes{1}(1)], 'color', shadeColor)
if numel(cyc.drugEndTimes{1}) > 1 
    plot_stim_shading([cyc.drugStartTimes{1}(2), cyc.drugEndTimes{1}(2)] - cyc.drugStartTimes{1}(1), ...
        'color', rgb('green'))
end
xlim([xx(1) xx(end)])
title('Vector strength')
% xlim([-90 300])

% Min - max
f = figure(2); clf; hold on;
f.Color = [1 1 1];
cm = lines(numel(unique(cyc.expID)));
allExpIDs = unique(cyc.expID);
xx = cycTimes - cyc.drugStartTimes{1}(1);
for iExp = 1:numel(allExpIDs)
    currExpID = allExpIDs{iExp};
    currCyc = cyc(strcmp(cyc.expID, currExpID), :);
    if singleRois
        for iRoi = 1:size(currCyc, 1)
            plot(xx, smoothdata(currCyc.cycRange{iRoi}, 1, 'gaussian', smWin), 'color', ...
                    cm(iExp, :));
%               jitRange = mean(diff(xx)) * 0.1;
%                 test = (2*jitRange) rand(numel(currCyc.cycVp{iRoi}), 1)
%               plot(xx, currCyc.cycVp{iRoi}, '*');
        end
    elseif singleExpts
        currData = smoothdata(cell2mat(currCyc.cycRange'), 1, 'gaussian', smWin);
        SE = std_err(currData, 2);
        plot(xx, mean(currData, 2), 'linewidth', 2, 'color', cm(iExp, :));
        legendStr{end + 1} = currCyc.expID{1};
        if shadeExpSEM
            upperY = (mean(currData, 2) + SE);
            lowerY = (mean(currData, 2) - SE);
            jbfill(xx', upperY', lowerY', cm(iExp, :), cm(iExp, :), 1, 0.2);
        end
    end
end
meanData = mean(smoothdata(cell2mat(cyc.cycRange'), 1, 'gaussian', smWin), 2);
plot(xx, meanData, 'color', 'k', 'linewidth', 3);
if numel(cyc.drugStartTimes{1}) > 1
    shadeColor = rgb('orange');
else
    shadeColor = rgb('green');
end
plot_stim_shading([0, cyc.drugEndTimes{1}(1) - cyc.drugStartTimes{1}(1)], 'color', shadeColor)
if numel(cyc.drugEndTimes{1}) > 1 
    plot_stim_shading([cyc.drugStartTimes{1}(2), cyc.drugEndTimes{1}(2)] - cyc.drugStartTimes{1}(1), ...
        'color', rgb('green'))
end
if shadeSEM
    SE = std_err(smoothdata(cell2mat(cyc.cycRange'), 1, 'gaussian', 1), 2);
    upperY = (meanData + SE);
    lowerY = (meanData - SE);
    jbfill(xx', upperY', lowerY', rgb('black'), rgb('black'), 1, 0.2);
end
title('Min - max dF/F')
xlim([xx(1) xx(end)])

% PVA offset SD
f = figure(3);clf;
hold on;
f.Color = [1 1 1];
cyc = allCycExpData;
cm = lines(numel(unique(cyc.expID)));
allExpIDs = unique(cyc.expID);
xx = cycTimes - cyc.drugStartTimes{1}(1);
legendStr = {};
if singleExpts
    for iExp = 1:numel(allExpIDs)
        currExpID = allExpIDs{iExp};
        currCyc = cyc(strcmp(cyc.expID, currExpID), :);
        currData = smoothdata(currCyc.cycOffsetStd{:}, 'gaussian', smWin);
        plot(xx, currData, 'linewidth', 2, 'color', cm(iExp, :));
    end
end
meanData = mean(smoothdata(cell2mat(cyc.cycOffsetStd'), 1, 'gaussian', smWin), 2);
plot(xx, meanData, 'color', 'k', 'linewidth', 3);
if shadeSEM
    SE = std_err(smoothdata(cell2mat(cyc.cycOffsetStd'), 1, 'gaussian', smWin), 2);
    upperY = (meanData + SE);
    lowerY = (meanData - SE);
    jbfill(xx', upperY', lowerY', rgb('black'), rgb('black'), 1, 0.2);
end
if numel(cyc.drugStartTimes{1}) > 1
    shadeColor = rgb('orange');
else
    shadeColor = rgb('green');
end
plot_stim_shading([0, cyc.drugEndTimes{1}(1) - cyc.drugStartTimes{1}(1)], 'color', shadeColor)
if numel(cyc.drugEndTimes{1}) > 1 
    plot_stim_shading([cyc.drugStartTimes{1}(2), cyc.drugEndTimes{1}(2)] - cyc.drugStartTimes{1}(1), ...
        'color', rgb('green'))
end
xlim([xx(1) xx(end)])
title('Bar-bump offset SD')

% Bump amplitude
f = figure(4);clf;
hold on;
f.Color = [1 1 1];
cm = lines(numel(unique(cyc.expID)));
allExpIDs = unique(cyc.expID);
xx = cycTimes - cyc.drugStartTimes{1}(1);
legendStr = {};
if singleExpts
    for iExp = 1:numel(allExpIDs)
        currExpID = allExpIDs{iExp};
        currCyc = cyc(strcmp(cyc.expID, currExpID), :);
        currData = smoothdata(currCyc.cycBumpAmp{:}, 'gaussian', smWin);
        plot(xx, currData, 'linewidth', 2, 'color', cm(iExp, :));
    end
end
meanData = mean(smoothdata(cell2mat(cyc.cycBumpAmp'), 1, 'gaussian', smWin), 2);
plot(xx, meanData, 'color', 'k', 'linewidth', 3);
if shadeSEM
    SE = std_err(smoothdata(cell2mat(cyc.cycBumpAmp'), 1, 'gaussian', smWin), 2);
    upperY = (meanData + SE);
    lowerY = (meanData - SE);
    jbfill(xx', upperY', lowerY', rgb('black'), rgb('black'), 1, 0.2);
end
if numel(cyc.drugStartTimes{1}) > 1
    shadeColor = rgb('orange');
else
    shadeColor = rgb('green');
end
plot_stim_shading([0, cyc.drugEndTimes{1}(1) - cyc.drugStartTimes{1}(1)], 'color', shadeColor)
if numel(cyc.drugEndTimes{1}) > 1 
    plot_stim_shading([cyc.drugStartTimes{1}(2), cyc.drugEndTimes{1}(2)] - cyc.drugStartTimes{1}(1), ...
        'color', rgb('green'))
end
xlim([xx(1) xx(end)])
title('Bump amplitude')













catch ME; rethrow(ME); end

%%
cyc = innerjoin(allCycRoiData, allCycExpData);
    
% epochs = [-90 0; ...
%            120, 180; ...
%            180, 240; ...
%            240 300]; ...
%       
epochs = [-90 0; ...
           60, 150; ...
           210 300];
% epochs = [-90, 0; ...
%            210,  300];   

epochs = epochs + 600;


% Get mean vector strength and phase within each epoch
epochVsMeans = [];
epochVpMeans = [];
for iRoi = 1:size(cyc, 1)
    t = cyc.cycTimes{iRoi};
    atpStart = cyc.drugStartTimes{iRoi}(1);
    for iEpoch = 1:size(epochs, 1)
        epochCycles = t > (epochs(iEpoch, 1) + atpStart) & t < (epochs(iEpoch, 2) + atpStart);
        epochVsMeans(iRoi, iEpoch) = mean(cyc.cycVs{iRoi}(epochCycles));
        epochVpMeans(iRoi, iEpoch) = circ_mean(cyc.cycVp{iRoi}(epochCycles));  
    end
end

% Get the mean bar-bump offset, offset SD, and bump amplitude for each epoch
epochOffsetMean = [];
epochOffsetSD = [];
epochBumpAmpMean = [];
for iExp = 1:numel(currExpList)
    currExpID = currExpList{iExp};
    currTrialNums = trialNums{iExp};
    currTbl = tbl(strcmp(tbl.expID, currExpID) & ismember(tbl.trialNum, currTrialNums), :);
    
    % Get volTimes for the whole block
    expVolTimes = [];
    totalTrialDur = 0;
    expDrugStartTimes = [];
    expDrugEndTimes = [];
    for iTrial = 1:size(currTbl, 1)
        currVolTimes = cell2mat(currTbl(currTbl.trialNum == currTrialNums(iTrial), :).volTimes);
        expVolTimes = [expVolTimes; currVolTimes + totalTrialDur];
        if ~isnan(currTbl.startTime(iTrial))
            expDrugStartTimes(end + 1) = currTbl.startTime(iTrial) + totalTrialDur;
            expDrugEndTimes(end + 1) = expDrugStartTimes(end) + currTbl.duration(iTrial);
        end
        totalTrialDur = totalTrialDur + currTbl(currTbl.trialNum == currTrialNums(iTrial), :).trialDuration;
    end
    
    % Extract bump offset, amplitude, and panels pos for entire exp (filling in volumes from 
    % darkness trials with NaN as needed)
    expOffset = []; expBumpAmp = []; expPanelsPos = [];
    for iTrial = 1:size(currTbl, 1)
        if currTbl.usingPanels(iTrial)
            expOffset = [expOffset; currTbl.pvaOffset{iTrial}];
            expBumpAmp = [expBumpAmp; currTbl.bumpAmp{iTrial}];
            expPanelsPos = [expPanelsPos; currTbl.panelsPosVols{iTrial}];
        else
            expOffset = [expOffset; nan(currTbl.nVolumes(iTrial), 1)];
            expBumpAmp = [expBumpAmp; nan(currTbl.nVolumes(iTrial), 1)];
            expPanelsPos = [expPanelsPos; nan(currTbl.nVolumes(iTrial), 1)];
        end
    end
    
    % Identify bar positions that are behind the fly and set them to NaN in the data
    barBehindFlyPanelsPositions = [1:8, 81:96] - 1; % Subtracting 1 for zero-indexed
%     barBehindFlyPanelsPositions = [1:12, 77:96] - 1; % Subtracting 1 for zero-indexing
    barBehindFlyVols = ismember(expPanelsPos, barBehindFlyPanelsPositions);
    expOffset(barBehindFlyVols) = nan;
    expBumpAmp(barBehindFlyVols) = nan;

    for iEpoch = 1:size(epochs, 1)
        relEpochTimes = epochs(iEpoch, :) + expDrugStartTimes(1);
        epochVols = expVolTimes > relEpochTimes(1) & expVolTimes < relEpochTimes(2);
        
        currOffset = expOffset(epochVols) + pi;
        currBumpAmp = expBumpAmp(epochVols);
        
        epochOffsetMean(iExp, iEpoch) = circ_mean(currOffset(~isnan(currOffset)));
        [~, epochOffsetSD(iExp, iEpoch)] = circ_std(currOffset(~isnan(currOffset)));
        epochBumpAmpMean(iExp, iEpoch) = mean(currBumpAmp, 'omitnan');
    end
end


% Bar-bump offset SD
f = figure(5); clf;
f.Color = [1 1 1];
hold on;
legendStr = {};
xx = 1:size(epochOffsetSD, 2);
for iExp = 1:size(epochOffsetSD, 1)
    plot(xx, epochOffsetSD(iExp, :), '-s', 'linewidth', 2);
    legendStr{end + 1} = currExpList{iExp};
end
legend(legendStr, 'autoupdate', 'off')
plot(xx, mean(epochOffsetSD, 1), '-s', 'color', 'k', ...
        'linewidth', 3, 'markersize', 12, 'markerfaceColor', 'k');
xlim([0.5 xx(end) + 0.5]);
title('Bar-bump offset SD')

% Mean bump offset
f = figure(6); clf;
f.Color = [1 1 1];
hold on;
legendStr = {};
xx = 1:size(epochOffsetMean, 2);
for iExp = 1:size(epochOffsetMean, 1)
    plot(xx, epochOffsetMean(iExp, :), '-s', 'linewidth', 2);
    legendStr{end + 1} = currExpList{iExp};
end
xlim([0.5 xx(end) + 0.5]);
legend(legendStr, 'autoupdate', 'off')
title('Mean bump offset')

% Mean delta offsets
f = figure(7); clf;
f.Color = [1 1 1];
hold on;
legendStr = {};
xx = 1:(size(epochOffsetMean, 2) - 1);
allYY = [];
for iExp = 1:size(epochOffsetMean, 1)
    for iEpoch = 1:(size(epochOffsetMean, 2) - 1)
        allYY(iExp, iEpoch) = abs(circ_dist(epochOffsetMean(iExp, iEpoch + 1), ...
                epochOffsetMean(iExp, iEpoch)));
    end
    plot(xx, allYY(iExp, :), '-s', 'linewidth', 2);
    legendStr{end + 1} = currExpList{iExp};
end
plot(xx, mean(allYY, 1), '-s', 'color', 'k', ...
        'linewidth', 3, 'markersize', 12, 'markerfaceColor', 'k');
xlim([0.5 xx(end) + 0.5]);
legend(legendStr, 'autoupdate', 'off')
title('Mean delta offset')
% ylim([0 2.5])

% Mean bump amplitude
f = figure(8); clf;
f.Color = [1 1 1];
hold on;
legendStr = {};
xx = 1:size(epochBumpAmpMean, 2);
for iExp = 1:size(epochBumpAmpMean, 1)
    plot(xx, epochBumpAmpMean(iExp, :), '-s', 'linewidth', 2);
    legendStr{end + 1} = currExpList{iExp};
end
plot(xx, mean(epochBumpAmpMean, 1), '-s', 'color', 'k', ...
        'linewidth', 3, 'markersize', 12, 'markerfaceColor', 'k');
xlim([0.5 xx(end) + 0.5]);
legend(legendStr, 'autoupdate', 'off')
title('Mean bump amplitude')



% relPhaseMeans = epochVpMeans + pi;
% % relPhaseMeans = relPhaseMeans - repmat(relPhaseMeans(:, 1), 1, size(relPhaseMeans, 2));
% 
% test = unwrap(relPhaseMeans, [], 2);
% yy = (test - repmat(test(:, 1), 1, size(test, 2)))';
% f = figure(7);clf; 
% violinplot(abs(diff(yy, 1, 1)'));
% title('relative vector phase shifts')
% 
% f = figure(8); clf;
% f.Color = [1 1 1];
% hold on;
% SEM = std_err(epochVsMeans, 1);
% for iRoi = 1:size(epochVsMeans, 1)
%     errorbar(1:size(epochVsMeans, 2), epochVsMeans(iRoi, :), epochVsSEM(iRoi, :), '-s');%, ...
% end
% errorbar(1:size(epochVsMeans, 2), mean(epochVsMeans, 1)', SEM, '-s', 'linewidth', 3, 'color', 'k', ...
%         'markerSize', 12, 'markerfacecolor', 'k')
% xlim([0.5 size(epochVsMeans, 2) + 0.5]);
% % ylim([0, max(relPhaseMeans(:))]);
% title('vector strength')











