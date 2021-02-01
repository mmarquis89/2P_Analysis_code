%% Extract single-cycle summary data for all experiments

% 5-5-5-5-5 timing
darkExpList = {'20201201-1', '20201203-1', '20201203-2', '20201210-1', '20201210-2', '20201201-2'};
darkExpTrialNums = {3:7, 3:7, 1:5, 3:7, 1:5, 2:6};

% % 5-5-10 timing
% visExpList = {'20201117-1', '20201117-3', '20201120-2', '20210118-1', '20210119-1'};
% visExpTrialNums = {1:4, 2:5, 2:5, 2:5, 2:6};

% % 5-5-5 timing
% visExpList = {'20201117-1', '20201117-3', '20201117-4', '20201120-2', '20201201-3', '20210118-1', ...
%     '20210118-2', '20210119-1'};
% visExpTrialNums = {1:3, 2:4, 2:4, 2:4, 3:5, 2:4, 2:4, 2:4};

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

currExpList = visExpList;
trialNums = visExpTrialNums;
% currExpList = darkExpList;
% trialNums = darkExpTrialNums;

allCycData = [];
for iExp = 1:numel(currExpList)
    
    currExpID = currExpList{iExp};
    currTrialNums = trialNums{iExp};   

    % Separate data for current experiment and trial(s)
    currTbl = tbl(strcmp(tbl.expID, currExpID) & ismember(tbl.trialNum, currTrialNums), :);
    currCycleData = cycleData(strcmp(cycleData.expID, currExpID) & ismember(cycleData.trialNum, ...
        currTrialNums), :);
    
    drugTrials = find(~isnan(currTbl.startTime));
    
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
                currTbl(drugTrials(iTrial), :).startTime;
            drugEndTimes(iTrial, 1) = drugStartTimes(iTrial) + currTbl(drugTrials(iTrial), ...
                    :).duration;
        end
    end
    
    % Create new row for each ROI and append to the appropriate table
    for iRoi = 1:size(cycVp, 2)
        newRow = table({currExpID}, iRoi, 'variableNames', {'expID', 'roiNum'});
        newRow.cycTimes = {cycTimes'};
        newRow.cycVs = {cycVs(:, iRoi)};
        newRow.cycVp = {cycVp(:, iRoi)};
        newRow.cycMin = {cycMins(:, iRoi)};
        newRow.cycMax = {cycMaxes(:, iRoi)};
        newRow.cycRange = {cycRanges(:, iRoi)};
        newRow.drugStartTimes = {drugStartTimes};
        newRow.drugEndTimes = {drugEndTimes};
        allCycData = [allCycData; newRow];
    end
    allCycData = [allCycData; newRow];
end

%

smWin = 3;
singleRois = 0;
singleExpts = 0;
shadeSEM = 1;
shadeExpSEM = 1;

cyc = allCycData;

cycTimes = cyc.cycTimes{1};
Vs = cell2mat(cyc.cycVs');
Vp = cell2mat(cyc.cycVp');
cycMin = cell2mat(cyc.cycMin');
cycMax = cell2mat(cyc.cycMax');
cycRanges = cell2mat(cyc.cycRange');

% Vector strength
figure(1); clf; hold on;
cm = lines(numel(unique(cyc.expID)));
allExpIDs = unique(cyc.expID);
xx = cycTimes - cyc.drugStartTimes{1}(1);
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
        if shadeExpSEM
            upperY = (mean(currData, 2) + SE);
            lowerY = (mean(currData, 2) - SE);
            jbfill(xx', upperY', lowerY', cm(iExp, :), cm(iExp, :), 1, 0.2);
        end
    end
end
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

% Min - max
figure(2); clf; hold on;
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
        end
    elseif singleExpts
        currData = smoothdata(cell2mat(currCyc.cycRange'), 1, 'gaussian', smWin);
        SE = std_err(currData, 2);
        plot(xx, mean(currData, 2), 'linewidth', 2, 'color', cm(iExp, :));
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
    SE = std_err(smoothdata(Vs, 1, 'gaussian', smWin), 2);
    upperY = (meanData + SE);
    lowerY = (meanData - SE);
    jbfill(xx', upperY', lowerY', rgb('black'), rgb('black'), 1, 0.2);
end
title('Min - max dF/F')
xlim([xx(1) xx(end)])















