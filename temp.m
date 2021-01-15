
%% Load data

% parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\new_PPL201_experiments';

expList = {'20201222-1', '20201222-2', '20201228-1', '20201228-2', '20201228-3', '20210102-1', ...
        '20210102-2', '20210102-3', '20210102-4'};

% expList = PPL201_expList(7:end);
    
% Load metadata 
[expMd, trialMd] = load_metadata(expList, parentDir);

% Load imaging data
roiData = load_roi_data(expList, parentDir);

% Load FicTrac data
ftData = load_ft_data(expList, parentDir);

% Load event data
eventData = load_event_data(expList, parentDir);
locEventData = eventData.locomotion;
% locEventData = eventData.flailing;
odorEventData = eventData.odor;

% Generate consolidated source data table
tbl = inner_join(roiData, expMd, trialMd);
tbl = outerjoin(tbl, ftData, 'type', 'left', 'mergekeys', 1);

sourceTbl = tbl;
%% 

roiName = 'TypeB';

trialNum = 2;

smWin = 5;

tbl = sourceTbl(contains(sourceTbl.roiName, roiName), :);
tbl = tbl(tbl.trialNum == trialNum, :);

tblExps = unique(tbl.expID);

for iExp = 1:numel(tblExps)
    
    currTbl = tbl(strcmp(tbl.expID, tblExps{iExp}), :);
    currTbl = currTbl(1, :);
    
    volTimes = (1:currTbl.nVolumes) ./ currTbl.volumeRate;
    currOdorEvents = odorEventData.eventData(strcmp(odorEventData.eventData.expID, tblExps{iExp}) & ...
            odorEventData.eventData.trialNum == trialNum, :);
    odorStimTimes = [currOdorEvents.onsetTime, currOdorEvents.offsetTime];
%     odorOnsetTimes = currTbl.optoStimTiming{:}(1):sum(currTbl.optoStimTiming{:}(2:3)):volTimes(end);
%     odorOffsetTimes = odorOnsetTimes + currTbl.optoStimTiming{:}(2);
%     odorStimTimes = [odorOnsetTimes; odorOffsetTimes]';
    
%     currTrialList = locEventData.trialList(strcmp(locEventData.trialList.expID, currTbl.expID), :);
    locVols = locEventData.create_logical_array(currTbl.nVolumes, volTimes, currTbl(:, {'expID', ...
            'trialNum'}));
    
    figure(iExp); clf; 
    
    ax = subaxis(2, 1, 1, 'ml', 0.05, 'mr', 0.05); hold on;
    yy = smoothdata(currTbl.rawFl{:}, 1, 'gaussian', smWin);
    yy(logical(locVols)) = nan;
    plot(volTimes, yy);
    plot_stim_shading(odorStimTimes);
    ax.XLim = [0 volTimes(end)];
    title(currTbl.expID{:})
    
    ax = subaxis(2, 1, 2); hold on;
    plot(currTbl.frameTimes{:}, repeat_smooth(currTbl.moveSpeed{:}, 10, 'smWin', smWin))
    plot_stim_shading(odorStimTimes);
    plot(volTimes, 10 * locVols, 'color', 'r')
    ax.XLim = [0 currTbl.frameTimes{:}(end)];
end


%% Group data from several compatible trials into a single block

expID = '20201222-1';

trialNums = [];

% Extract data from current experiment
currExpData = tbl(strcmp(tbl.expID, expID), :);
if ~isempty(trialNums)
   currExpData = currExpData(ismember(currExpData.trialNum, trialNums), :); 
end
nTrials = size(unique(currExpData.trialNum), 1);
roiNames = unique(currExpData.roiName);
nRois = numel(roiNames);

% Put all fl data in a single matrix for each fl type
rawFlMat = zeros(currExpData.nVolumes(1), nTrials, nRois);
trialDffMat = rawFlMat;
expDffMat = rawFlMat;
for iRoi = 1:nRois
    currRoiData = currExpData(strcmp(currExpData.roiName, roiNames{iRoi}), :); 
    for iTrial = 1:nTrials
        rawFlMat(:, iTrial, iRoi) = currRoiData.rawFl{iTrial}; % [volume, trial, ROI]
        trialDffMat(:, iTrial, iRoi) = (currRoiData.rawFl{iTrial} - ...
                currRoiData.trialBaseline(iTrial)) ./ currRoiData.trialBaseline(iTrial); % [volume, trial, ROI]
        expDffMat(:, iTrial, iRoi) = (currRoiData.rawFl{iTrial} - currRoiData.expBaseline(1)) ./ ...
                currRoiData.expBaseline(1); % [volume, trial, ROI]
    end
end
volTimes = (1:size(rawFlMat, 1)) ./ currExpData.volumeRate(1);
nVolumes = numel(volTimes);

%%

%% PLOT SINGLE-TRIAL HEATMAPS FOR ENTIRE BLOCK

saveFig = 0;

smWin = 5;
flType = 'expDff';
roiName = 'ANT';

figPos = [1800 950];

% Adjust plot spacing and margins
SV = 0.002;
SH = 0.02;
ML = 0.04;
MR = 0.02;
MT = 0.05;
MB = 0.065;



% Get correct source data
targetRoiInd = find(contains(roiNames, roiName));
if strcmp(flType, 'rawFl')
    flData = rawFlMat(:, :, targetRoiInd);
elseif strcmp(flType, 'trialDff')
    flData = trialDffMat(:, :, targetRoiInd);
elseif strcmp(flType, 'expDff')
    flData = expDffMat(:, :, targetRoiInd);
end
plotFl = smoothdata(flData, 1, 'gaussian', smWin, 'omitnan'); % --> [volume, trial]

% To give all figures the same color scale
plotFl(1, 1, :) = min(plotFl(:), [], 'omitnan');
plotFl(end, end, :) = max(plotFl(:), [], 'omitnan');

f = figure(3);clf;
f.Color = [1 1 1];
if ~isempty(figPos)
    f.Position(3:4) = figPos;
    if sum(f.Position([2, 4])) > 1080
        f.Position(2) = 50;
    end
end
colormap(magma);
for iTrial = 1:nTrials 
    
    % Plot Fl Data
    subaxis(nTrials, 8, [1:7] + (8 * (iTrial - 1)), 'mt', MT, 'mb', MB, 'sv', SV, 'mr', MR, ...
            'ml', ML, 'sh', SH);
    imagesc([0, currRoiData.trialDuration(1)], [1, 2], plotFl(:, iTrial)');
    hold on
    
    % Label each plot with trial number
    ax = gca();
    ax.YTickLabel = [];
    t = ylabel(['Trial #', num2str(currRoiData.trialNum(iTrial))], 'rotation', 0, 'FontSize', 12);
    t.HorizontalAlignment = 'right';
    t.VerticalAlignment = 'middle';
    t.Position(1) = t.Position(1) * 2;
    
    % Label X axis on final trial
    if iTrial < nTrials
        ax.XTickLabel = [];
    else
        str = ax.YLabel.String;
        pos = ax.YLabel.Position;
        ax.FontSize = 12;
        xlabel('Time (sec)', 'FontSize', 14);
        ylabel(str, 'rotation', 0, 'FontSize', 12);      
    end
    
end%iTrial















