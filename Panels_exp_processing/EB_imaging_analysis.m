
%% Set directories and get expList
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\EB-DAN_GroupedAnalysisData';
figDir = fullfile(parentDir, 'Figs');

fileList = dir(fullfile(parentDir, '2020*'));
expList = unique(cellfun(@(x) x(1:10), {fileList.name}, 'uniformoutput', 0));

%% Load data for the selected experiments

% Load metadata 
[expMd, trialMd] = load_metadata(expList, parentDir);

% Load imaging data
roiData = load_roi_data(expList, parentDir);

% Load FicTrac data
ftData = load_ft_data(expList, parentDir);

% Load flailing event data
try
    eventData = load_event_data(expList, parentDir);
    flailingEvents = eventData.flailing;
catch; end

% Load panels metadata
panelsMetadata = load_panels_metadata(expList, parentDir);

%% GENERATE MOVE SPEED ANALYSIS OBJECT 

a = MoveSpeedAnalysis(expList', parentDir);

%% RUN MOVE SPEED ANALYSIS

saveFig = 1;

expNums = [];

skipTrials = {[11:15], [11:16], [8:15], [1:4, 9, 12:13], 7:10, 8:12}';

useSavedObj = 0;

speedType = 'yaw';

plotting_minSpeed = 0;
nAnalysisBins = 40;

roiName = 'FB-DAN';
flType = 'expDff'; % 'slidingbaseline', 'normalized', 'trialDff', 'expDff', 'rawFl' 
figDims = [1520 840];

try
    
xLim = [0 360];
yLim = [];
fontSize = 12;

% Adjust plot spacing and margins
SV = 0.08;
SH = 0.03;
ML = 0.05;
MR = 0.02;
MT = 0.04;
MB = 0.07;

if isempty(expNums)
   expNums = 1:numel(expList);
else
    skipTrials = skipTrials(expNums);
end
if numel(expNums) == 3
    plotLayoutDims = [2, 2];
else
    plotLayoutDims = numSubplots(numel(expNums));
end
savedAnalysisObjects = [];
% if isempty(savedAnalysisObjects) || ~useSavedObj
%     savedAnalysisObjects = allPlotParams;
%     savedAnalysisObjects.obj = repmat({[]}, size(allPlotParams, 1), 1);
% end

f = figure(1); clf;
f.Color = [1 1 1];
if ~isempty(figDims)
   f.Position(3:4) = figDims; 
end
allRs = [];
for iExp = 1:numel(expNums)   
    
    % Get current expID and params
    currExpID = expList{expNums(iExp)};
    disp(currExpID);
    
    a.params.roiName = roiName;
    a.params.nHistBins = 25;
    a.params.flType = flType;
    a.params.convertSpeedUnits = 0;
    a.params.smWinVols = 5;
    a.params.smWinFrames = 7;
    a.params.nSmoothRepsFrames = 20;
    a.params.skipTrials = skipTrials{iExp};
    a.params.plotting_minSpeed = plotting_minSpeed;
    a.params.nAnalysisBins = nAnalysisBins;
    
    % Create axes for current plot
    ax = subaxis(plotLayoutDims(1), plotLayoutDims(2), iExp, ...
            'ml', ML, 'mr', MR, 'mt', MT, 'mb', MB, 'sv', SV, 'sh', SH);
    
    if ~useSavedObj || isempty(savedAnalysisObjects.obj{expNums(iExp)})
        
        % Initialize filters
        a.filterDefs.expID = currExpID;
        a.filterDefs.roiName = roiName;
        a = a.init_filters();
        
        % If the filter matches more than one ROI, just use the first one
        uniqueRois = unique(a.sourceDataTable.subset.roiName);
        if numel(uniqueRois) > 1
            a.filterDefs.roiName = uniqueRois{1};
            a = a.init_filters();
        end
        
        % Run analysis
        a = a.analyze();
        
        savedAnalysisObjects.obj{expNums(iExp)} = a;
    else
        % Re-use existing analysis object
        a = savedAnalysisObjects.obj{expNums(iExp)};
    end
%     
%     % Plot binned speed vs. dF/F
%     a = a.generate_binned_flData(flType, speedType);
%     a.plot_binned_fl(ax);
%     ax.FontSize = fontSize;
%     ax.YLabel.String = 'Mean dF/F';
%     if strcmp(speedType, 'forward')
%         cf = a.analysisOutput.r(1);
%     else
%         cf = a.analysisOutput.r(2);
%     end
%     allRs(end + 1) = cf;
%     ax.Title.String = [currExpID, ' — ', roiName, ' — r = ', num2str(cf, 2)];
%     box off
%     yStr = flType;

%     % Forward speed vs. Yaw
%     [binMidpoints, binMeans, binSEM] = MoveSpeedAnalysis.bin_data(a.analysisOutput.fwSpeedData, ...
%             a.analysisOutput.yawSpeedData, a.params.nAnalysisBins, a.params.plotting_minSpeed);
%     errorbar(ax, binMidpoints, binMeans, binSEM, '-o', 'color', 'k', 'linewidth', 1);
%     ax.XLabel.String = 'Forward speed (mm/sec)';
%     ax.YLabel.String = 'Mean yaw speed (deg/sec)';
%     ax.FontSize = fontSize;
%     ax.Title.String = [currExpID, ' — ', roiName];    
%     yStr = 'yawSpeed';
%     box off
% 
     % Yaw vs. forward speed
    [binMidpoints, binMeans, binSEM] = MoveSpeedAnalysis.bin_data(a.analysisOutput.yawSpeedData, ...
            a.analysisOutput.fwSpeedData, a.params.nAnalysisBins, a.params.plotting_minSpeed);
    errorbar(ax, binMidpoints, binMeans, binSEM, '-o', 'color', 'k', 'linewidth', 1);
    ax.XLabel.String = 'Yaw speed (deg/sec)';
    ax.YLabel.String = 'Mean forward speed (mm/sec)';
    ax.FontSize = fontSize;
    ax.Title.String = [currExpID, ' — ', roiName];   
    yStr = 'fwSpeed';
    box off

    % Remove axis labels if not in last row (for X) or first column (for Y)
    if iExp <= numel(expNums) - plotLayoutDims(2)
        ax.XLabel.String = [];
    end
    if ~ismember(iExp, 1:plotLayoutDims(2):numel(expNums))
        ax.YLabel.String = [];
    end

    % Manually set axis limits if necessary
    if ~isempty(xLim)
        xlim(ax, xLim);
    end
    if ~isempty(yLim)
        ylim(ax, yLim);
    end
end%iExp

if saveFig
    if isempty(xLim)
        xlStr = 'scaledX';
    else
        xlStr = 'matchedX';
    end
    if isempty(yLim)
        ylStr = 'scaledY';
    else
        ylStr = 'matchedY';
    end
    fileName = [roiName, '_binned_', speedType, 'Speed_vs_', yStr, '_', xlStr, '_', ylStr]
    save_figure(f, figDir, fileName);
end


catch ME; rethrow(ME); end

%% Plot ROI data

smWin = 5;
smReps = 50;

% expNums = [];
% trialNums = repmat({[]}, numel(expList), 1);
expNums = [7];
trialNums = {[], []};
% trialNums = {[1:6, 11:16], [1:7, 13:17]};

% roiNames, moveSpeed, fwSpeed, yawSpeed
ax1_varNames = {'FB-DAN', 'fwSpeed'};
ax2_varNames = {'fwSpeed', 'yawSpeed'};


allVarNames = [ax1_varNames, ax2_varNames];
roiNames = unique(allVarNames(ismember(allVarNames, unique(roiData.roiName))));

cm = custom_colormap(numel(ax1_varNames) + numel(ax2_varNames));
cm = [0, 0, 0; cm];

if isempty(expNums)
   expNums = 1:numel(expList); 
end
currExpList = expList;
currExpList = currExpList(expNums);
for iExp = 1:numel(currExpList)
        
    % Get current ExpID and trial num
    currExpID = currExpList{iExp};
    if isempty(trialNums{iExp})
        currTrialNums = 1:max(roiData(strcmp(roiData.expID, currExpID), :).trialNum);
    else
        currTrialNums = trialNums{iExp};
    end
    currTrialMd = trialMd(strcmp(trialMd.expID, currExpID), :);
    
    % Extract ROI data
    currRoiData = roiData(strcmp(roiData.expID, currExpID), :);
    flMat = []; flData = [];
    for iROI = 1:numel(roiNames)
        currFl = cell2mat(currRoiData(strcmp(currRoiData.roiName, roiNames{iROI}), :).rawFl');
        flMat = currFl(:, currTrialNums);
        flData(:, iROI) = as_vector(smoothdata(flMat, 1, 'gaussian', smWin));
    end

    % Extract FicTrac data
    allMoveSpeed = []; allFwSpeed = []; allYawSpeed = [];
    for iTrial = 1:numel(currTrialNums)
        currFtData = ftData(strcmp(ftData.expID, currExpID) & ...
                ftData.trialNum == currTrialNums(iTrial), :);
        moveSpeed = repeat_smooth(currFtData.moveSpeed{:}, smReps, 'smWin', smWin);
        fwSpeed = repeat_smooth(currFtData.fwSpeed{:}, smReps, 'smWin', smWin);
        yawSpeed = abs(repeat_smooth(currFtData.yawSpeed{:}, smReps, 'smWin', smWin));
        sideSpeed = repeat_smooth(currFtData.sideSpeed{:}, smReps, 'smWin', smWin);
        totalSpeed = (abs(fwSpeed) + abs(yawSpeed) + abs(sideSpeed)) ./ 3;
        frameRate = 1/ median(diff(currFtData.frameTimes{:}));
        
        % Get downsampled FicTrac data
        volTimes = (1:currTrialMd(currTrialMd.trialNum == currTrialNums(iTrial), :).nVolumes) ...
                ./ expMd(strcmp(expMd.expID, currExpID), :).volumeRate;
        moveSpeedVols = []; fwSpeedVols = []; yawSpeedVols = [];
        for iVol = 1:numel(volTimes)
            dsFrame = argmin(abs(currFtData.frameTimes{:} - volTimes(iVol)));
            moveSpeedVols(iVol) = moveSpeed(dsFrame);
            fwSpeedVols(iVol) = fwSpeed(dsFrame);
            yawSpeedVols(iVol) = yawSpeed(dsFrame);
        end
        
        % Get panels X position data
        nPanelsFrames = double(currTrialMd.nPanelsFrames(iTrial));
        panelsFrameTimes = (1:nPanelsFrames) ./ ...
                expMd(strcmp(expMd.expID, currExpID), :).panelsDisplayRate;
        panelsPosX = panelsMetadata
        
        allMoveSpeed = [allMoveSpeed; smoothdata(moveSpeedVols', 1, 'gaussian', smWin)];
        allFwSpeed = [allFwSpeed; smoothdata(fwSpeedVols', 1, 'gaussian', smWin)];
        allYawSpeed = [allYawSpeed; smoothdata(yawSpeedVols', 1, 'gaussian', smWin)];
        
    end%iTrial
    
    % Create source data table, normalizing each variable so its max == 1
    tbl = table(allMoveSpeed ./ max(allMoveSpeed), allFwSpeed ./ max(allFwSpeed), allYawSpeed ...
            ./ max(allYawSpeed), 'variablenames', {'moveSpeed', 'fwSpeed', 'yawSpeed'});
    for iRoi = 1:numel(roiNames)
        currFl = flData(:, iRoi) - min(flData(:, iRoi));
        tbl.(roiNames{iRoi}) = currFl ./ max(currFl);
    end
    
    % Create figure
    f = figure(iExp);clf; hold on;
    f.Color = [1 1 1];
    
    % Plot 1
    if isempty(ax2_varNames)
        nPlots = 1;
    else
        nPlots = 2;
    end
    ax1 = subaxis(nPlots, 1, 1, 'sv', 0.04, 'm', 0.05);
    hold on;
    cmCount = 1;
    xx = (1:size(flData, 1)) ./ expMd(strcmp(expMd.expID, currExpID), :).volumeRate;
    for iVar = 1:numel(ax1_varNames)
        plot(xx, tbl.(ax1_varNames{iVar}), 'linewidth', 1.5, 'color', cm(cmCount, :));
        cmCount = cmCount + 1;
    end
    xlim([xx(1), xx(end)])
    legend(ax1_varNames);
    
    % Plot 2
    if nPlots == 2
        ax2 = subaxis(2, 1, 2);
        hold on;
        xx = (1:size(flData, 1)) ./ expMd(strcmp(expMd.expID, currExpID), :).volumeRate;
        for iVar = 1:numel(ax2_varNames)
            plot(xx, tbl.(ax2_varNames{iVar}), 'linewidth', 1.5, 'color', cm(cmCount, :));
            cmCount = cmCount + 1;
        end
        xlim([xx(1), xx(end)])
        legend(ax2_varNames);
    end
    
    % Link X-axis limits
    linkaxes([ax1, ax2], 'x')    
    
    
%     % yyaxis right
%     plot(xx, fwSpeedVols ./ max(fwSpeedVols), 'linewidth', 2, 'color', 'b')
%     hold on;
%     plotYaw = smoothdata(abs(yawSpeedVols), 2, 'gaussian', smWin);
%     plot(xx, plotYaw ./ max(plotYaw), '-', 'linewidth', 2, 'color', 'r')
%     legend([roiNames, {'fwSpeed', 'yawSpeed'}])

    % % ROI 1 + ROI 2
    % allFl = cell2mat(roiData(strcmp(roiData.roiName, roiNames{1}), :).rawFl');
    % ax = subaxis(3, 1, 1, 'sv', 0.04, 'm', 0.05);
    % plot((1:numel(allFl(:, trialNums)))/numel(allFl(:, trialNums)), smoothdata(allFl(:, trialNums), 1, 'gaussian', smWin), 'linewidth', 1);
    % hold on;
    % % yyaxis right
    % allFl = cell2mat(roiData(strcmp(roiData.roiName, roiNames{2}), :).rawFl');
    % plot((1:numel(allFl(:, trialNums)))/numel(allFl(:, trialNums)), smoothdata(allFl(:, trialNums), 1, 'gaussian', smWin), 'linewidth', 1);
    % 
    % % ROI 3
    % allFl = cell2mat(roiData(strcmp(roiData.roiName, roiNames{3}), :).rawFl');
    % ax = subaxis(3, 1, 3);
    % plot(smoothdata(allFl(:, trialNums), 1, 'gaussian', smWin), 'linewidth', 1);
    % 
    % % % Panels
    % % ax = subaxis(3, 1, 2);
    % % % allPosX = cell2mat(panelsMetadata.panelsPosX)
    % % plot(1:numel(panelsMetadata.panelsPosX{1}), as_vector(panelsMetadata.panelsPosX{1}'), ...
    % %         'linewidth', 1)
    % %   
% 
%     % Plot 2
%     ax2 = subaxis(2, 1, 2);
%     xx = (1:numel(moveSpeed))' ./ frameRate;
%     plot(xx, fwSpeed, 'linewidth', 2); hold on;
%     xlim([xx(1), xx(end)])
%     ylim([0 20])
%     yyaxis right
%     plot(xx, abs(rad2deg(yawSpeed)), 'linewidth', 2);
%     legend({'Fw', 'Yaw'})
%     % yyaxis right
%     % plot((1:numel(moveSpeed))/numel(moveSpeed), abs(yawSpeed), 'linewidth', 1);
%     % plot((1:numel(totalSpeed))/numel(totalSpeed), totalSpeed, 'linewidth', 1);
%     % plot((1:numel(fwSpeed))/numel(fwSpeed), abs(fwSpeed), 'linewidth', 1);
%     linkaxes([ax1, ax2], 'x')
    %
% 
%     % Calculate fl-speed correlations for each ROI
%     for iRoi = 1:size(flMat, 2)
%         moveR = corrcoef(flMat(:, iRoi), moveSpeedVols');
%         fwR = corrcoef(flMat(:, iRoi), abs(fwSpeedVols)');
%         yawR = corrcoef(flMat(:, iRoi), abs(yawSpeedVols)');
% 
%         disp('  ')
%         disp(roiNames{iRoi})
%         disp(['move: ', num2str(moveR(2,1), 2)])
%         disp(['fw: ', num2str(fwR(2,1), 2)])
%         disp(['yaw: ', num2str(yawR(2,1), 2)])
%     end

end

%% Pairwise scatter plots of all ROIs

currRoiData = roiData(strcmp(roiData.expID, '20201027-1'), :);

try
roiNames = unique(currRoiData.roiName);

roiNames = roiNames();

nROIs = numel(roiNames);
f = figure(1);clf; hold on;
f.Color = [1 1 1];
ax = subaxis(nROIs, nROIs, 1, 'sv', 0.04, 'm', 0.05);
count = 1;
for i = 1:nROIs
    ROI_1 = roiNames{i};
    disp(ROI_1)
    fl_1 = cell2mat(currRoiData(strcmp(currRoiData.roiName, ROI_1), :).rawFl);
    for j = 1:nROIs
        ROI_2 = roiNames{j};
        fl_2 = cell2mat(currRoiData(strcmp(currRoiData.roiName, ROI_2), :).rawFl);
        if i ~= j 
            subaxis(nROIs, nROIs, count)
            scatter(fl_1, fl_2);
            title([ROI_1, ' vs ' ROI_2])
        end
        count = count + 1;
    end
end
catch ME; rethrow(ME); end
