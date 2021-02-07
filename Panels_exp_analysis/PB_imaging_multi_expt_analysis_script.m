%% Extract single-cycle summary data for all experiments

% 5-5-5-5-5 timing
darkExpList = {'20201201-1', '20201203-1', '20201203-2', '20201210-1', '20201210-2', '20201201-2'};
darkExpTrialNums = {3:7, 3:7, 1:5, 3:7, 1:5, 2:6};

% % 5-5-10 timing
% visExpList = {'20201117-1', '20201117-3', '20201120-2', '20210118-1', '20210119-1'};
% visExpTrialNums = {1:4, 2:5, 2:5, 2:5, 2:6};

% % 5-5-5 timing
% visExpList = {'20201117-1', '20201117-3', '20201117-4', '20201120-2', '20210118-1', ...
%     '20210118-2', '20210119-1', '20201201-3'};
% visExpTrialNums = {1:3, 2:4, 2:4, 2:4, 2:4, 2:4, 2:4, 3:5};

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
% % 
% 5-5-5 vis genetic control
visExpList = {'20201120-1', '20201120-3'};
visExpTrialNums = {2:4, 2:4};

visExpTrialNums = visExpTrialNums(~strcmp(visExpList, '20210118-2'));
visExpList = visExpList(~strcmp(visExpList, '20210118-2'));

currExpList = visExpList;
trialNums = visExpTrialNums;
% currExpList = darkExpList;
% trialNums = darkExpTrialNums;

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

saveFig = 0;
fileNameSuffix = '_withVisualStim';

smWin = 3;
singleRois = 0;
singleExpts = 1;
shadeSEM = 0;
shadeExpSEM = 0;


% For visual stim experiments
epochs = [-120 0; 120 240; 360 480];% -2-0 min, 2-4 min, 6-8 min
% 
% % For darkness+ATP experiments
% lateEpochs = epochs + 600;
% epochs = [epochs(1, :); 300 420; lateEpochs];
% % epochs(1:2, :) = epochs(1:2, :) + 105; 
% epochs(1, :) = epochs(1, :) - 45; % Because there's no data after the ATP trial starts at t=(-45)


disp(epochs ./ 60)

try
cyc = innerjoin(allCycRoiData, allCycExpData);

cycTimes = cyc.cycTimes{1};
Vs = cell2mat(cyc.cycVs');
Vp = cell2mat(cyc.cycVp');
cycMin = cell2mat(cyc.cycMin');
cycMax = cell2mat(cyc.cycMax');
cycRanges = cell2mat(cyc.cycRange');
% 
% % Vector strength
% f = figure(1); clf; hold on;
% f.Color = [1 1 1];
% cm = lines(numel(unique(cyc.expID)));
% allExpIDs = unique(cyc.expID);
% xx = cycTimes - cyc.drugStartTimes{1}(1);
% legendStr = {};
% for iExp = 1:numel(allExpIDs)
%     currExpID = allExpIDs{iExp};
%     currCyc = cyc(strcmp(cyc.expID, currExpID), :);
%     if singleRois
%         for iRoi = 1:size(currCyc, 1)
%             yy = smoothdata(currCyc.cycVs{iRoi}, 1, 'gaussian', smWin);
%             plot(xx(~isnan(yy)), yy(~isnan(yy)) , 'color', ...
%                     cm(iExp, :));
%         end
%     elseif singleExpts
%         currData = smoothdata(cell2mat(currCyc.cycVs'), 1, 'gaussian', smWin);
%         SE = std_err(currData, 2);
%         plot(xx, mean(currData, 2), 'linewidth', 2, 'color', cm(iExp, :));
%         legendStr{end + 1} = currCyc.expID{1};
%     end
% end
% % legend(legendStr, 'autoupdate', 'off');
% meanData = mean(smoothdata(Vs, 1, 'gaussian', smWin), 2);
% plot(xx, meanData, 'color', 'k', 'linewidth', 3);
% if shadeSEM
%     SE = std_err(smoothdata(Vs, 1, 'gaussian', smWin), 2);
%     upperY = (meanData + SE);
%     lowerY = (meanData - SE);
%     jbfill(xx', upperY', lowerY', rgb('black'), rgb('black'), 1, 0.2);
% end
% if numel(cyc.drugStartTimes{1}) > 1
%     shadeColor = rgb('orange');
% else
%     shadeColor = rgb('green');
% end
% plot_stim_shading([0, cyc.drugEndTimes{1}(1) - cyc.drugStartTimes{1}(1)], 'color', shadeColor)
% if numel(cyc.drugEndTimes{1}) > 1 
%     plot_stim_shading([cyc.drugStartTimes{1}(2), cyc.drugEndTimes{1}(2)] - cyc.drugStartTimes{1}(1), ...
%         'color', rgb('green'))
% end
% xlim([xx(1) xx(end)])
% title('Vector strength')
% yL = ylim();
% if ~isempty(epochs)
%     for iEpoch = 1:size(epochs, 1)
%         plot(epochs(iEpoch, :), [1 1] * yL(2) * 0.98, 'linewidth', 8, 'color', 'k')
%     end
% end
% 
% % Min - max
% f = figure(2); clf; hold on;
% f.Color = [1 1 1];
% cm = lines(numel(unique(cyc.expID)));
% allExpIDs = unique(cyc.expID);
% xx = cycTimes - cyc.drugStartTimes{1}(1);
% for iExp = 1:numel(allExpIDs)
%     currExpID = allExpIDs{iExp};
%     currCyc = cyc(strcmp(cyc.expID, currExpID), :);
%     if singleRois
%         for iRoi = 1:size(currCyc, 1)
%             plot(xx, smoothdata(currCyc.cycRange{iRoi}, 1, 'gaussian', smWin), 'color', ...
%                     cm(iExp, :));
%         end
%     elseif singleExpts
%         currData = smoothdata(cell2mat(currCyc.cycRange'), 1, 'gaussian', smWin);
%         SE = std_err(currData, 2);
%         plot(xx, mean(currData, 2), 'linewidth', 2, 'color', cm(iExp, :));
%         legendStr{end + 1} = currCyc.expID{1};
%     end
% end
% meanData = mean(smoothdata(cell2mat(cyc.cycRange'), 1, 'gaussian', smWin), 2);
% plot(xx, meanData, 'color', 'k', 'linewidth', 3);
% if numel(cyc.drugStartTimes{1}) > 1
%     shadeColor = rgb('orange');
% else
%     shadeColor = rgb('green');
% end
% plot_stim_shading([0, cyc.drugEndTimes{1}(1) - cyc.drugStartTimes{1}(1)], 'color', shadeColor)
% if numel(cyc.drugEndTimes{1}) > 1 
%     plot_stim_shading([cyc.drugStartTimes{1}(2), cyc.drugEndTimes{1}(2)] - cyc.drugStartTimes{1}(1), ...
%         'color', rgb('green'))
% end
% if shadeSEM
%     SE = std_err(smoothdata(cell2mat(cyc.cycRange'), 1, 'gaussian', 1), 2);
%     upperY = (meanData + SE);
%     lowerY = (meanData - SE);
%     jbfill(xx', upperY', lowerY', rgb('black'), rgb('black'), 1, 0.2);
% end
% title('Min - max dF/F')
% xlim([xx(1) xx(end)])
% yL = ylim();
% if ~isempty(epochs)
%     for iEpoch = 1:size(epochs, 1)
%         plot(epochs(iEpoch, :), [1 1] * yL(2) * 0.98, 'linewidth', 8, 'color', 'k')
%     end
% end


% PVA offset SD
f = figure(3);clf;
hold on;
f.Color = [1 1 1];
f.Position(3:4) = [945 495];
cyc = allCycExpData;
cm = lines(numel(unique(cyc.expID)));
allExpIDs = unique(cyc.expID);
xx = cycTimes - cyc.drugStartTimes{1}(1);
legendStr = {};
if singleExpts
    for iExp = 1:numel(allExpIDs)
        currExpID = allExpIDs{iExp};
        legendStr{end + 1} = currExpID;
        currCyc = cyc(strcmp(cyc.expID, currExpID), :);
        currData = rad2deg(smoothdata(currCyc.cycOffsetStd{:}, 'gaussian', smWin));
%         plot(xx, currData, 'linewidth', 1, 'color', cm(iExp, :));
        plot(xx(~isnan(currData)), currData(~isnan(currData)), 'linewidth', 1, 'color', 0.5*[1 1 1])
    end
end
meanData = rad2deg(mean(smoothdata(cell2mat(cyc.cycOffsetStd'), 1, 'gaussian', smWin), 2));
plot(xx, meanData, 'color', 'k', 'linewidth', 5);
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
ylim([0 180]);
yL = ylim();
% plot_stim_shading([0, cyc.drugEndTimes{1}(1) - cyc.drugStartTimes{1}(1)], 'color', shadeColor)
plot([0, cyc.drugEndTimes{1}(1) - cyc.drugStartTimes{1}(1)], [1 1] * yL(2) * 0.98, 'color', ...
        shadeColor, 'linewidth', 8)
if numel(cyc.drugEndTimes{1}) > 1 
%     plot_stim_shading([cyc.drugStartTimes{1}(2), cyc.drugEndTimes{1}(2)] - cyc.drugStartTimes{1}(1), ...
%         'color', rgb('green'))
    plot([cyc.drugStartTimes{1}(2), cyc.drugEndTimes{1}(2)] - cyc.drugStartTimes{1}(1), ...
            [1 1] * yL(2) * 0.98, 'color', rgb('green'), 'linewidth', 8)
end
xlim([xx(1) xx(end)])
ylabel('Bar-bump offset SD')
% legend(legendStr, 'autoupdate', 'off')
if ~isempty(epochs)
    for iEpoch = 1:size(epochs, 1)
        plot(epochs(iEpoch, :), [1 1] * yL(2) * 0.98, 'linewidth', 8, 'color', 'k')
    end
end
ax = gca();
ax.FontSize = 14;
xlabel('Time (sec)');
if saveFig
    f.UserData.expList = currExpList;
    f.UserData.trialNums = trialNums;
    f.UserData.smWin = smWin;
    fileName = ['EB-DAN_+ATP_mean_bar-bump_offset_SD', fileNameSuffix];
    save_figure(f, figDir, fileName);
end



% 
% % Bump amplitude
% f = figure(4);clf;
% hold on;
% f.Color = [1 1 1];
% cm = lines(numel(unique(cyc.expID)));
% allExpIDs = unique(cyc.expID);
% xx = cycTimes - cyc.drugStartTimes{1}(1);
% legendStr = {};
% if singleExpts
%     for iExp = 1:numel(allExpIDs)
%         currExpID = allExpIDs{iExp};
%         currCyc = cyc(strcmp(cyc.expID, currExpID), :);
%         currData = smoothdata(currCyc.cycBumpAmp{:}, 'gaussian', smWin);
%         plot(xx, currData, 'linewidth', 2, 'color', cm(iExp, :));
%     end
% end
% meanData = mean(smoothdata(cell2mat(cyc.cycBumpAmp'), 1, 'gaussian', smWin), 2);
% plot(xx, meanData, 'color', 'k', 'linewidth', 3);
% if shadeSEM
%     SE = std_err(smoothdata(cell2mat(cyc.cycBumpAmp'), 1, 'gaussian', smWin), 2);
%     upperY = (meanData + SE);
%     lowerY = (meanData - SE);
%     jbfill(xx', upperY', lowerY', rgb('black'), rgb('black'), 1, 0.2);
% end
% if numel(cyc.drugStartTimes{1}) > 1
%     shadeColor = rgb('orange');
% else
%     shadeColor = rgb('green');
% end
% plot_stim_shading([0, cyc.drugEndTimes{1}(1) - cyc.drugStartTimes{1}(1)], 'color', shadeColor)
% if numel(cyc.drugEndTimes{1}) > 1 
%     plot_stim_shading([cyc.drugStartTimes{1}(2), cyc.drugEndTimes{1}(2)] - cyc.drugStartTimes{1}(1), ...
%         'color', rgb('green'))
% end
% xlim([xx(1) xx(end)])
% title('Bump amplitude')
% yL = ylim();
% if ~isempty(epochs)
%     for iEpoch = 1:size(epochs, 1)
%         plot(epochs(iEpoch, :), [1 1] * yL(2) * 0.98, 'linewidth', 8, 'color', 'k')
%     end
% end

catch ME; rethrow(ME); end

% Plot mean values across experimental epochs
cyc = innerjoin(allCycRoiData, allCycExpData);

% Get mean vector strength and phase within each epoch
epochVsMeans = [];
epochVpMeans = [];
for iRoi = 1:size(cyc, 1)
    t = cyc.cycTimes{iRoi};
    atpStart = cyc.drugStartTimes{iRoi}(1);
    for iEpoch = 1:size(epochs, 1)
        epochCycles = t > (epochs(iEpoch, 1) + atpStart) & t < (epochs(iEpoch, 2) + atpStart);
        epochVsMeans(iRoi, iEpoch) = mean(cyc.cycVs{iRoi}(epochCycles));
        epochVpMeans(iRoi, iEpoch) = circ_mean(cyc.cycVp{iRoi}(epochCycles) + pi);  
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
xx = 1:size(epochOffsetSD, 2)-1;
for iExp = 1:size(epochOffsetSD, 1)
    plot(xx, rad2deg(epochOffsetSD(iExp, 1:end-1)), '-', 'linewidth', 2);
    legendStr{end + 1} = currExpList{iExp};
end
% legend(legendStr, 'autoupdate', 'off')
plot(xx, rad2deg(mean(epochOffsetSD(:, 1:end-1), 1)), '-s', 'color', 'k', ...
        'linewidth', 3, 'markersize', 12, 'markerfaceColor', 'k');
xlim([0.5 xx(end) + 0.5]);
ylabel('Bar-bump offset SD')
ax = gca();
ax.XTick = xx;
ax.FontSize = 12;
ax.XTickLabel = {'Before', 'After'};
title('ATP + visual', 'fontsize', 12)
ylim([0 120])
f.Position(3:4) = [220 400];
if saveFig
    f.UserData.smWin = smWin;
    f.UserData.epochs = epochs;
    f.UserData.expList = currExpList;
    f.UserData.trialNums = trialNums;
    fileName = ['EB-DAN_prePost_ATP_offsetSD_crunch', fileNameSuffix];
    save_figure(f, figDir, fileName);
end

% % Mean bump offset
% f = figure(6); clf;
% f.Color = [1 1 1];
% hold on;
% legendStr = {};
% xx = 1:size(epochOffsetMean, 2);
% for iExp = 1:size(epochOffsetMean, 1)
%     plot(xx, epochOffsetMean(iExp, :), '-', 'linewidth', 2);
%     legendStr{end + 1} = currExpList{iExp};
% end
% xlim([0.5 xx(end) + 0.5]);
% % legend(legendStr, 'autoupdate', 'off')
% title('Mean bump offset')
% ax = gca();
% ax.XTick = xx;
% ylim([-pi, pi]);

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
    plot(xx, rad2deg(allYY(iExp, :)), '-', 'linewidth', 2);
    legendStr{end + 1} = currExpList{iExp};
end
plot(xx, rad2deg(mean(allYY, 1)), '-s', 'color', 'k', ...
        'linewidth', 3, 'markersize', 12, 'markerfaceColor', 'k');
xlim([0.5 xx(end) + 0.5]);
% legend(legendStr, 'autoupdate', 'off')
ylabel('Mean \Delta offset')
ax = gca();
ax.XTick = xx;
ylim([0 120]);
ax.FontSize = 12;
title('ATP + visual', 'fontsize', 12)
ax.XTickLabel = {'+ATP', 'post-ATP'};
f.Position(3:4) = [220 400];
if saveFig
    f.UserData.smWin = smWin;
    f.UserData.epochs = epochs;
    f.UserData.expList = currExpList;
    f.UserData.trialNums = trialNums;
    fileName = ['EB-DAN_prePost_ATP_meanDeltaOffset_crunch', fileNameSuffix];
    save_figure(f, figDir, fileName);
end

% % Mean bump amplitude
% f = figure(8); clf;
% f.Color = [1 1 1];
% hold on;
% legendStr = {};
% xx = 1:size(epochBumpAmpMean, 2);
% for iExp = 1:size(epochBumpAmpMean, 1)
%     plot(xx, epochBumpAmpMean(iExp, :), '-', 'linewidth', 2);
%     legendStr{end + 1} = currExpList{iExp};
% end
% plot(xx, mean(epochBumpAmpMean, 1), '-s', 'color', 'k', ...
%         'linewidth', 3, 'markersize', 12, 'markerfaceColor', 'k');
% xlim([0.5 xx(end) + 0.5]);
% % legend(legendStr, 'autoupdate', 'off')
% title('Mean bump amplitude')
% ax = gca();
% ax.XTick = xx;
% 
% % Vector strength
% f = figure(9); clf;
% f.Color = [1 1 1];
% hold on;
% xx = 1:size(epochVsMeans, 2);
% plot(xx, epochVsMeans, '-', 'linewidth', 2);
% plot(xx, mean(epochVsMeans, 1), '-s', 'color', 'k', ...
%         'linewidth', 3, 'markersize', 12, 'markerfaceColor', 'k');
% xlim([0.5 xx(end) + 0.5]);
% title('Vector strength');
% ax = gca();
% ax.XTick = xx;

% Mean preferred bar position shifts for individual EB wedges
f = figure(10); clf;
f.Color = [1 1 1];
hold on;
xx = 1:(size(epochVpMeans, 2) - 1);
allYY = [];
for iRoi = 1:size(epochVpMeans, 1)
    for iEpoch = 1:(size(epochVpMeans, 2) - 1)
        allYY(iRoi, iEpoch) = abs(circ_dist(epochVpMeans(iRoi, iEpoch + 1), ...
                epochVpMeans(iRoi, iEpoch)));
    end
    plot(xx, rad2deg(allYY(iRoi, :)), '-', 'linewidth', 2);
end
% violinplot(allYY, []);
plot(xx, rad2deg(mean(allYY, 1)), '-s', 'color', 'k', ...
        'linewidth', 3, 'markersize', 12, 'markerfaceColor', 'k');
xlim([0.5 xx(end) + 0.5]);
ylabel('Mean \Delta preferred bar position')
ax = gca();
ax.XTick = xx;
ylim([0 180]);
ax.FontSize = 12;
title('ATP + visual', 'fontsize', 12)
ax.XTickLabel = {'+ATP', 'post-ATP'};
f.Position(3:4) = [220 400];
if saveFig
    f.UserData.smWin = smWin;
    f.UserData.epochs = epochs;
    f.UserData.expList = currExpList;
    f.UserData.trialNums = trialNums;
    fileName = ['EB-DAN_prePost_ATP_meanDeltaPreferredPos_crunch', fileNameSuffix];
    save_figure(f, figDir, fileName);
end




