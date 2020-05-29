
%% CREATE AND SAVE EventAlignedData OBJECTS

saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Saved_AlignEvent_objects';

% Create a temporary EventAlignedData object just to get a list of the event types 
alignObj = EventAlignedData(load_expList());

% Get list of all available event types
eventNames = fieldnames(alignObj.eventObjects);
eventNames = eventNames(~strcmp(eventNames, 'quiescence')); % Prob don't want to analyze this

% Process and save an object for each event type
for iType = 2:numel(eventNames)
    
    currEventType = eventNames{iType};
    
    % Reload alignObj if it's already been used
    if iType > 1
        alignObj = EventAlignedData(load_expList());
    end
    
    % Set alignment event, then load and process data
    alignObj = alignObj.set_align_event(currEventType);
    alignObj = alignObj.extract_roi_data();
    alignObj = alignObj.extract_fictrac_data();
    alignObj = alignObj.create_filter_event_vectors();
    
    % Save object
    saveDate = datestr(datetime(), 'yyyymmdd');
    saveFileName = fullfile(saveDir, [saveDate, '_AlignEventObj_', currEventType, '.mat']);
    save(saveFileName, 'alignObj');
    
end

%%

test3 = alignObj.output_analysis_subset([2 4]);
unique(test3.roiName)

%%

filts = struct();
filts.roiName = 'Background';

filts.ballstop = 0;
filts.flailing = 0;
filts.grooming = 0;
filts.isolatedmovement = 1;
filts.locomotion = 0;
filts.odor = 0;
filts.optostim = 0;
filts.panelsflash = 0;
filts.quiescence = nan;
filts.soundstim = 0;

% filts.odorName = 'EtOH';
% filts.concentration = '4%'

filtNames = fieldnames(filts);
filterVec = ones(size(test3, 1), 1);
for iFilt = 1:size(filtNames)
    currFiltName = filtNames{iFilt};
    if ~isnan(filts.(currFiltName))
        if isa(filts.(currFiltName), 'char')
            % Filter based on a string field
            currFiltVec = ~cellfun(@isempty, regexp(test3.(currFiltName), filts.(currFiltName), ...
                'once'));
        elseif strcmp(currFiltName, test3.Properties.CustomProperties.alignEventName)
            % Filter out events with pre-onset instances of the alignment event type
            onsetVols = cellfun(@(x) argmin(abs(x)), test3.volTimes);
            currFiltVec = ~cellfun(@(x, y) any(x(1:(y - 1))), test3.(currFiltName), ...
                    mat2cell(onsetVols, ones(size(onsetVols))));
            
        elseif filts.(currFiltName) == 0
            % Filter out based on an event filter vector
            currFiltVec = ~cellfun(@any, test3.(currFiltName));
        end
        filterVec = filterVec .* currFiltVec;
    end
end
filterVec = logical(filterVec);
test4 = test3(filterVec, :);

%

test4 = innerjoin(test4, oneNoteMd(:, {'expID', 'olfactometerVersion'}));

smWin = 3;
singleTrials = 1;

expIDList = unique(test4.expID);
subplotDims = numSubplots(numel(expIDList));
if all(subplotDims == [1 3])
   subplotDims = [2 2]; 
end
f = figure(1); clf;
f.Color = [1 1 1];
currPlotNum = 0;
for iExp = 1:numel(expIDList)
    currExpID = expIDList{iExp};
    currData = test4(strcmp(test4.expID, currExpID), :);
    
    % If more than one ROI matched the filter string, only use the first one
    firstRoiName = currData.roiName{1};
    currData = currData(strcmp(currData.roiName, firstRoiName), :);
    
    currEventFl = cell2mat(currData.eventFl');
    
    fl = currEventFl;
    bl = repmat(currData.trialBaseline', size(fl, 1), 1);
    currEventFl = (fl - bl) ./ bl;
    
    currEventFl = currEventFl(:, ~any(isnan(currEventFl), 1));
    
    [~, onsetVol] = min(abs(currData.volTimes{1}));
    stimDur = currData.offsetTime(1) - currData.onsetTime(1);
    [~, offsetVol] = min(abs(currData.volTimes{1} - stimDur));

%     currEventFl((onsetVol - 1):(offsetVol + 1), :) = nan;
%     currEventFl((onsetVol - 1):(onsetVol + 2), :) = nan;

    %    currEventMoveSpeed = cell2mat(currData.moveSpeed');
    if size(currEventFl, 2) > 0
        
        currPlotNum = currPlotNum + 1;
        ax = subaxis(subplotDims(1), subplotDims(2), currPlotNum, 'm', 0.02, 'sv', 0.04, 'sh', 0.02);
        if singleTrials
            xx = repmat(currData.volTimes{1}', 1, size(currEventFl, 2));
            plot(xx, smoothdata(currEventFl, 1, 'gaussian', smWin, 'omitnan'));
        end
        hold on;
        plot(currData.volTimes{1}, mean(smoothdata(currEventFl, 1, 'gaussian', smWin), 2, 'omitnan'), ...
                'linewidth', 3, 'color', 'k');

        meanData = mean(smoothdata(currEventFl, 1, 'gaussian', smWin), 2, 'omitnan');
        SE = std_err(currEventFl, 2);
        xx = currData.volTimes{1};
        upperY = (meanData + SE)';
        lowerY = (meanData - SE)';
        jbfill(xx(~isnan(upperY)), upperY(~isnan(upperY)), lowerY(~isnan(lowerY)), [0 0 0], [0 0 0], 1, 0.2);
        yL = ylim();
        
%         plot_stim_shading([onsetVol, offsetVol]);
        plot([0, 0], yL, 'color', 'r', 'linewidth', 3);
        
% imagesc(isnan(currEventFl'));
        ylim(yL);
%                imagesc(currEventFl');
        % imagesc(currEventMoveSpeed')
        title([currExpID, ' - Ver', num2str(currData.olfactometerVersion(1))]);
        box off
    end
end

unique(test3.roiName)
% unique(test3(:, {'odorName', 'concentration'}))


%%
singleExpData = test4(strcmp(test4.expID, '20171108-1'), :);
plotVars = {'eventFl', 'moveSpeed', 'locomotion', 'isolatedmovement', 'panelsflash', 'ballstop'};


plotVars = plotVars([1 4 6]);
removeNanTrials = 1;


for iVar = 1:numel(plotVars)
    
    dataArr = cell2mat(singleExpData.(plotVars{iVar})');
    if removeNanTrials
       dataArr = dataArr(:, ~any(isnan(dataArr), 1));
    end
    
    figure(100 + iVar); clf; 
    imagesc(dataArr');
    title(plotVars{iVar})
end


%%


%
% 
% 
% 
% 










