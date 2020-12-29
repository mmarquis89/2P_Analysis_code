
%% CREATE AND SAVE EventAlignedData OBJECTS

saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\new_PPL201_experiments\alignEvent_objects';
expList = {'20201222-1', '20201222-2'};

% Create a temporary EventAlignedData object just to get a list of the event types 
alignObj = EventAlignedData(expList, 'dataDir', 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\new_PPL201_experiments');

% Get list of all available event types
eventNames = fieldnames(alignObj.eventObjects);
eventNames = eventNames(~strcmp(eventNames, 'quiescence')); % Prob don't want to analyze this
 
% Process and save an object for each event type
for iType = 1:numel(eventNames)
    
    currEventType = eventNames{iType};
    disp(currEventType);
    
    % Reload alignObj if it's already been used
    if iType > 1
        alignObj = EventAlignedData(expList, 'dataDir', 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\new_PPL201_experiments');
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
