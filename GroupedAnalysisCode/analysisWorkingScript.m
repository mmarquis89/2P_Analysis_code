% % 
% %% Load a base aligned data object
% 
% parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Saved_AlignEvent_objects';
% alignEventDateStr = '20200528';
% 
% % ballstop, flailing, grooming, isolatedmovement, locomotion, odor, optostim, panelsflash, soundstim
% alignEventName = 'isolatedmovement';
% 
% load(fullfile(parentDir, [alignEventDateStr, '_AlignEventObj_', alignEventName, '.mat']));
% 
% %% Generate aligned event data table
% 
% analysisWin = [2 3];
% 
% dt = alignObj.output_analysis_DataTable(analysisWin);
% 
% %% Generate filter structure and use it to get a subset of the data
% filterDefs = alignObj.create_filterDefs();
% 
% % General fields
% filterDefs.expID = [];
% filterDefs.trialNum = [];
% filterDefs.expName = [];
% filterDefs.moveSpeed = [];
% filterDefs.yawSpeed = [];
% filterDefs.roiName = [];
% 
% % Event filter vectors
% filterDefs.ballstop = 0;
% filterDefs.grooming = 0;
% filterDefs.isolatedmovement = -1;
% filterDefs.locomotion = 0;
% filterDefs.flailing = 0;
% filterDefs.odor = 0;
% filterDefs.optostim = 0;
% filterDefs.panelsflash = 0;
% filterDefs.soundstim = 0;
% 
% % OneNote experiment metadata
% filterDefs.genotype = [];
% filterDefs.ageHrs = [];
% filterDefs.foodType = [];
% filterDefs.bedtime = [];
% filterDefs.starvationTimeHrs = [];
% filterDefs.olfactometerVersion = [];
% filterDefs.bathTemp = [];
% 
% % Alignment event-specific fields
% if strcmp(alignObj.alignEventName, 'odor')
%     filterDefs.odorName = [];
%     filterDefs.concentration = [];
% end
% % 
% dt = dt.generate_filters(filterDefs);
% % locTest = cell2mat(dt.subset.locomotion');
% % flTest = cell2mat(dt.subset.eventFl');
% % moveSpeedTest = cell2mat(dt.subset.moveSpeed');
% % figure; imagesc(locTest')
% % figure; imagesc(flTest')
% % figure; imagesc(moveSpeedTest')


%% PLOT FL DATA FOR EXPERIMENT, COLORED BY ROI

smWin = 1;
singleTrials = 0;
shadeStim = 0;
includeNaN = 0;
shadeSEM = 1;
minTrials = 5;

overlayVar = 'roiName';
overlayVarList = {'TypeD', 'TypeB', 'TypeF'};
overlayVarColors = [rgb('blue'); rgb('red'); rgb('green')];

transectionExpts = {'20190211-3', '20190216-1', '20190218-1', '20190218-2', '20190219-1', ...
        '20190219-2', '20190220-1', '20190226-1', '20190226-2', '20190226-3'};

overlayVarFilter = '';
for i=1:numel(overlayVarList)
   overlayVarFilter = [overlayVarFilter, '|', overlayVarList{i}];
end
overlayVarFilter = overlayVarFilter(2:end);
filterDefs.(overlayVar) = overlayVarFilter;
dt = dt.generate_filters(filterDefs);
expIDList = unique(dt.subset.expID);
subplotDims = numSubplots(numel(expIDList));
f = figure(1); clf;
f.Color = [1 1 1];
currPlotNum = 0;
clear ax;
emptyAxes = zeros(size(expIDList));
for iExp = 1:numel(expIDList)
    
    currPlotNum = currPlotNum + 1;
    ax(iExp) = subaxis(subplotDims(1), subplotDims(2), currPlotNum, 'm', 0.02, 'sv', 0.05, ...
            'sh', 0.03);
    hold on;
    emptyAxes(iExp) = 1;
    
    currExpID = expIDList{iExp};
    if ~ismember(currExpID, transectionExpts)
    dt = dt.update_filter('expID', currExpID);
%     dt = dt.update_filter(overlayVar, filterDefs.(filterIterVar));
%     iterVarList = unique(dt.subset.(filterIterVar));
    
    for iVar = 1:numel(overlayVarList)
        
        dt = dt.update_filter(overlayVar, overlayVarList{iVar});
        
        if size(dt.subset, 1) >= minTrials
            
            disp(unique(dt.subset.roiName)');
            
            if ~isempty(dt.subset)
                volTimes = dt.subset.volTimes{1};
            end
            
            fl = cell2mat(dt.subset.eventFl');
            bl = repmat(mean(fl(volTimes < 0, :), 1, 'omitnan'), size(fl, 1), 1);
%             bl = repmat(dt.subset.trialBaseline', size(fl, 1), 1);
            dff = (fl - bl) ./ bl;
            
            if ~includeNaN
                dff = dff(:, ~any(isnan(dff), 1));
            end
            
            % Get stim timing info
            [~, onsetVol] = min(abs(volTimes));
            stimDur = dt.subset.offsetTime(1) - dt.subset.onsetTime(1);
            [~, offsetVol] = min(abs(volTimes - stimDur));
            
%                    dff((onsetVol - 1):(offsetVol + 4), :) = nan;
%                    dff((onsetVol - 1):(onsetVol + 2), :) = nan;
            
%            Plot individual trials if appropriate
            if singleTrials
                xx = repmat(volTimes', 1, size(dff, 2));
                plot(xx, smoothdata(dff, 1, 'gaussian', smWin, 'omitnan'), 'color', overlayVarColors(iVar, :));
            end
            
            % Plot mean response for current ROI
            meanData = mean(smoothdata(dff, 1, 'gaussian', smWin), 2, 'omitnan');
            plot(volTimes', meanData, 'linewidth', 2, 'color', overlayVarColors(iVar, :));
            
            % Shade SEM if appropriate
            SE = std_err(dff, 2);
            xx = volTimes;
            upperY = (meanData + SE)';
            lowerY = (meanData - SE)';
            jbfill(xx(~isnan(upperY)), upperY(~isnan(upperY)), lowerY(~isnan(lowerY)), ...
                overlayVarColors(iVar, :), ...
                overlayVarColors(iVar, :), 1, 0.2);
            
            if any(~isnan(meanData)) && size(dff, 2) >= minTrials
                emptyAxes(iExp) = 0;
            end
            
        end%if
    end%iVar
            
    if ~emptyAxes(iExp)
        
        % Plot line or shade to mark stim onset and/or offset
        yL = ylim(); xL = xlim();
        if shadeStim
            plot_stim_shading([0 stimDur]);
        else
            plot([0, 0], yL, 'color', 'r', 'linewidth', 3);
        end
        ylim(yL); xlim(xL);
        
        title([currExpID]);
        box off
    end
    %         legend(roiList);

    % imagesc(isnan(currEventFl'));
    %                imagesc(currEventFl');
    % imagesc(currEventMoveSpeed')
    end
end%iExp

% Copy the non-empty plots over to a new figure
newFig = figure(2);
newFig.Color = [1 1 1];
clf;
subplotDims = numSubplots(sum(~emptyAxes));
if all(subplotDims == [1 3])
   subplotDims = [2 2]; 
end
goodAxes = ax(~emptyAxes);
for iAx = 1:numel(goodAxes)
    newAx = subaxis(subplotDims(1), subplotDims(2), iAx, 'ml', 0.04 , 'mt', 0.04, 'mb', 0.06, ...
            'mr', 0.02, 'sv', 0.05, 'sh', 0.03);
    axis off
    newPlot = copyobj(goodAxes(iAx), newFig);
    newPlot.Position = newAx.Position;
    newPlot.FontSize = 11;
    
    if mod(iAx, subplotDims(2)) == 1
        newPlot.YLabel.String = 'dF/F';
    end
    if iAx > (subplotDims(2) * (subplotDims(1) - 1))
        newPlot.XLabel.String = 'Time (sec)';
    end
end
% close(f);

% %% Save current figure
% 
% saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\testFigs';
% fileName = 'odorOnset_EtOH_neat_response_noMove_TypeD_TypeB_TypeF';
% 
% save_figure(newFig, saveDir, fileName)

















