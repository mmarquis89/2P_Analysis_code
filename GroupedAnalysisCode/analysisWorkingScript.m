%% Load a base aligned data object

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\';
alignEventDateStr = '20200609';

% ballstop, flailing, grooming, isolatedmovement, locomotion, odor, optostim, panelsflash, soundstim
alignEventName = 'odor';

load(fullfile(parentDir, 'Saved_AlignEvent_objects', [alignEventDateStr, '_AlignEventObj_', ...
        alignEventName, '.mat']));

%% Generate or load aligned event data table

analysisWin = [7 5];

saveFileName = [alignObj.alignEventName, '_pre_', num2str(analysisWin(1)), '_post_', ...
        num2str(analysisWin(2)), '.mat'];
if exist(fullfile(parentDir, 'Saved_DataTable_objects', saveFileName), 'file')
    disp('Loading saved DataTable...')
    load(fullfile(parentDir, 'Saved_DataTable_objects', saveFileName), 'dt');
    disp('DataTable loaded')
else
    disp('Generating new DataTable...')
    dt = alignObj.output_analysis_DataTable(analysisWin);
end

%% Generate filter structure and use it to get a subset of the data
baseFilterDefs = alignObj.create_filterDefs();

% General fields
baseFilterDefs.expID = '20180623-1';%'20180414-1|20180623-3|20181120-1|2019031503';
baseFilterDefs.trialNum = [];
baseFilterDefs.expName = [];
baseFilterDefs.moveSpeed = [];
baseFilterDefs.yawSpeed = [];
baseFilterDefs.roiName = 'TypeD';

% Event filter vectors
baseFilterDefs.ballstop = 0;
baseFilterDefs.grooming = [];
baseFilterDefs.isolatedmovement = [];
baseFilterDefs.locomotion = 0;
baseFilterDefs.flailing = 0;
baseFilterDefs.odor = -1;
baseFilterDefs.optostim = 0;
baseFilterDefs.panelsflash = 0;
baseFilterDefs.soundstim = 0;

% OneNote experiment metadata
baseFilterDefs.genotype = [];
baseFilterDefs.ageHrs = [];
baseFilterDefs.foodType = [];
baseFilterDefs.bedtime = [];
baseFilterDefs.starvationTimeHrs = [];
baseFilterDefs.olfactometerVersion = [];
baseFilterDefs.bathTemp = [];

% Alignment event-specific fields
if strcmp(alignObj.alignEventName, 'odor')
    baseFilterDefs.odorName = ['^.(?!.*background)(?!.*Oil)(?!.*Stop)(?!.*Air)'];%[];%
    baseFilterDefs.concentration = [];
end
% 
dt = dt.initialize_filters(baseFilterDefs);
baseFilterMat = dt.filterMat;


% SETUP PLOTTING VARIABLES

% Will make one subplot per value in this variable list
plotVar = 'expID';
transectionExpts = {'20190211-3', '20190216-1', '20190218-1', '20190218-2', '20190219-1', ...
        '20190219-2', '20190220-1', '20190226-1', '20190226-2', '20190226-3'};
expIDList = unique(dt.subset.expID);
plotVarList = expIDList(~ismember(expIDList, transectionExpts));

% expIDsOfInterest = {'20180207-2', '20180207-3', '20180209-1', '20180228-1', '20180228-3', '20180307-3', ...
%         '20180314-1', '20180316-2', '20180323-1', ...
%         '20180328-2', '20180329-1', '20180426-1', '20180525-1', '20180623-1', ...
%         '20180623-2', '20180627-1'};
% plotVarList = expIDsOfInterest';

% plotVar = 'roiName';
% plotVarList = {'VLP-AMMC'; 'TypeF'};



% Will plot one overlay group per value in this variable list (if rows matching the filter exist)

% groupVar = 'expID';
% transectionExpts = {'20190211-3', '20190216-1', '20190218-1', '20190218-2', '20190219-1', ...
%         '20190219-2', '20190220-1', '20190226-1', '20190226-2', '20190226-3'};
% expIDList = unique(dt.subset.expID)';
% groupVarList = expIDList(~ismember(expIDList, transectionExpts));
% groupVarColors = lines(numel(groupVarList));

% groupVar = 'roiName';
% groupVarList =  unique(dt.subset.roiName)';%{'TypeD', 'ANT'}; % {'TypeF', 'VLP-AMMC'};%{'TypeD', 'TypeB', 'TypeF'}; %
groupVar = 'odorName';
groupVarList = unique(dt.subset.odorName)';
groupVarColors = custom_colormap(numel(groupVarList));% [rgb('blue'); rgb('red'); rgb('green')]; %[rgb('green'); rgb('magenta')]; %



% Pre-calculate filter vectors so it's not being done repeatedly in a loop during plotting  
plotVarFilters = {};
groupVarFilters = {};

% Plot var
disp('Creating plotVar filters...')
for iPlot = 1:numel(plotVarList)
    plotVarFilters{iPlot} = dt.generate_filter(plotVar, plotVarList{iPlot});
end

% Grouping overlay var
disp('Creating grouping var filters...')
for iGroup = 1:numel(groupVarList)
    groupVarFilters{iGroup} = dt.generate_filter(groupVar, groupVarList{iGroup});
end

plotVarFilterTable = table(plotVarList, plotVarFilters', 'VariableNames', {plotVar, ...
    'filterVec'});
groupVarFilterTable = table(groupVarList', groupVarFilters', 'VariableNames', ...
        {groupVar, 'filterVec'});
    
disp('Plot and group filters ready')    

%% GENERATE PLOTS

% disp(unique(dt.subset.odorName))

% Set plotting options
smWin = 5;
singleTrials = 0;
shadeStim = 1;
includeNaN = 0;
shadeSEM = 1;
minTrials = 1;
matchYLims = 0;
useLegend = 1;
manualYLims = [];
minPlotGroups = 0;

% Configure plot spacing and margins
SV = 0.05;
SH = 0.025;
ML = 0.04;
MR = 0.015;
MT = 0.03;
MB = 0.06;

try 
    
% Get plot and group numbers
nPlots = size(plotVarFilterTable, 1);
nGroups = size(groupVarFilterTable, 1);

% Initialize figure 
f = figure(1); clf;
f.Color = [1 1 1];
clear ax;
subplotDims = numSubplots(nPlots);
emptyAxes = zeros(nPlots, 1);
globalYLims = [];
nPlottedGroups = zeros(nPlots, 1);
for iPlot = 1:nPlots
 
    ax(iPlot) = subaxis(subplotDims(1), subplotDims(2), iPlot);
    hold on;
    emptyAxes(iPlot) = 1;
    
    allStimDurs = [];
    for iGroup = 1:nGroups
        
        % Reset filter mat and filterDefs 
        dt.filterMat = baseFilterMat;
        dt.filterDefs = baseFilterDefs;
        
        % Manually set current filter status
        dt = dt.update_filter(plotVar, plotVarFilterTable.filterVec{iPlot});
        dt = dt.update_filter(groupVar, groupVarFilterTable.filterVec{iGroup});
        
        currSubset = dt.subset;
        
        if size(currSubset, 1) >= minTrials
            
            if ~isempty(currSubset)
                volTimes = currSubset.volTimes{1};
            end
            nPlottedGroups(iPlot) = nPlottedGroups(iPlot) + 1;
            
            fl = cell2mat(currSubset.eventFl');
            bl = repmat(mean(fl(volTimes < 0, :), 1, 'omitnan'), size(fl, 1), 1);
%             bl = repmat(dt.subset.trialBaseline', size(fl, 1), 1);
            dff = (fl - bl) ./ bl;
            
            if ~includeNaN
                dff = dff(:, ~any(isnan(dff), 1));
            end
            
            % Get stim timing info
            [~, onsetVol] = min(abs(volTimes));
            allStimDurs = [allStimDurs; currSubset.offsetTime - currSubset.onsetTime];
            
%                    dff((onsetVol - 1):(onsetVol + 2), :) = nan;
            
            % Plot individual trials if appropriate
            if singleTrials
                xx = repmat(volTimes', 1, size(dff, 2));
                plot(xx, smoothdata(dff, 1, 'gaussian', smWin, 'omitnan'), 'color', groupVarColors(iGroup, :));
            end
            
            % Plot mean response for current group
            meanData = mean(smoothdata(dff, 1, 'gaussian', smWin), 2, 'omitnan');
            plot(volTimes', meanData, 'linewidth', 2, 'color', groupVarColors(iGroup, :), ...
                    'displayname', groupVarList{iGroup});
            
            % Shade SEM if appropriate
            if shadeSEM
                SE = std_err(dff, 2);
                xx = volTimes;
                upperY = (meanData + SE)';
                lowerY = (meanData - SE)';
                jbfill(xx(~isnan(upperY)), upperY(~isnan(upperY)), lowerY(~isnan(lowerY)), ...
                    groupVarColors(iGroup, :), ...
                    groupVarColors(iGroup, :), 1, 0.2);
            end
            
            if any(~isnan(meanData)) && size(dff, 2) >= minTrials
                emptyAxes(iPlot) = 0;
                titleStr = [currSubset.expID{1}];
            end
            
        end%if
    end%iGroup
    
    if ~emptyAxes(iPlot)
        
        % Plot line or shade to mark stim onset and/or offset
        yL = ylim(); xL = xlim();
        if isempty(globalYLims)
            globalYLims = yL;
        else
            globalYLims(1) = min([globalYLims(1), yL(1)]);
            globalYLims(2) = max([globalYLims(2), yL(2)]);
        end
        if isempty(shadeStim)
            % Automatically decide whether to shade or just draw a line 
            if numel(unique(round(allStimDurs))) == 1
                shadeStim = 1;
            else
                shadeStim = 0;
            end
        end    
        if shadeStim
            plot_stim_shading([0 unique(allStimDurs(end))]);
        else
            plot([0, 0], [-100, 100], 'color', 'r', 'linewidth', 3); % Huge in case I change yLims later
        end
        ylim(yL); 
        xlim(xL);
        title([titleStr]);
        box off
    end

end%iExp

% Copy the non-empty plots over to a new figure
newFig = figure(2);
newFig.Color = [1 1 1];
clf;

emptyAxes(nPlottedGroups < minPlotGroups) = 1;

subplotDims = numSubplots(sum(~emptyAxes));
subplotDims(subplotDims == 0) = 1;
if all(subplotDims == [1 3])
   subplotDims = [2 2]; 
end
goodAxes = ax(~emptyAxes);
clear newAxes
for iAx = 1:numel(goodAxes)

    tempAx = subaxis(subplotDims(1), subplotDims(2), iAx, 'ml', ML, 'mr', MR, 'mt', MT, 'mb', MB, ...
            'sv', SV, 'sh', SH);
    axis off
    newAxes(iAx) = copyobj(goodAxes(iAx), newFig);
    newAxes(iAx).Position = tempAx.Position;
    newAxes(iAx).FontSize = 11;
    
    if mod(iAx, subplotDims(2)) == 1
        newAxes(iAx).YLabel.String = 'dF/F';
    end
    if iAx > (subplotDims(2) * (subplotDims(1) - 1))
        newAxes(iAx).XLabel.String = 'Time (sec)';
    end
    if ~isempty(manualYLims)
        ylim(newAxes(iAx), manualYLims);
    elseif matchYLims
        ylim(newAxes(iAx), globalYLims);
    end
    
   if iAx == numel(goodAxes) && useLegend
        handles = legendUnq(newFig);
        legend(newAxes(end), handles, 'fontsize', 14, 'location', 'best', 'autoupdate', 'off');
    end    
end
% close(f);

% % Plot legend on an empty axes in an unused part of the figure 
% tempAx = subaxis(subplotDims(1), subplotDims(2), numel(goodAxes) + 1, 'ml', ML, 'mr', MR, 'mt', MT, 'mb', MB, ...
%             'sv', SV, 'sh', SH);
% if  useLegend
%     handles = legendUnq(newFig);
%     legend(tempAx, handles, 'location', 'best', 'autoupdate', 'off');
%     set(gca, 'FontSize', 14)
%     axis off
% end

% Reset filter mat and filterDefs
dt.filterMat = baseFilterMat;
dt.filterDefs = baseFilterDefs;

catch ME; rethrow(ME); end

%% Save current figure

saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Figs';
fileName = 'groomingOnset_response_TypeD';

save_figure(newFig, saveDir, fileName)

















