function eventObjects = load_event_data(expList, parentDir)
% ==================================================================================================   
%  Loads event data files of ALL types for a set of expIDs and returns the event objects collected
%  in a struct. 
%  
%  INPUTS: 
%       expList                 = cell array of expIDs (YYYYMMDD-expNum format) to load data for.
%                                 Can provide a table with a field named "expID" instead.
%
%       parentDir (optional)    = override default location for source files
%
%  OUTPUTS:
%       eventObjects            = struct with each field containing an event object of the type
%                                 indicated by the field name, populated with all event dat of that 
%                                 type for the expIDs in expList
%
% ==================================================================================================

% Parse/process input arguments
if nargin < 2
    parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';
end
if isa(expList, 'table')
    expList = expList.expID;
end

% Identify all event data files for the selected experiments
eventDataFiles = [];
for iExp = 1:numel(expList)
    currExpID = expList{iExp};
    eventFiles = dir(fullfile(parentDir, [currExpID, '_event_data*.csv']));
    if ~isempty(eventFiles)
        eventDataFiles = [eventDataFiles; eventFiles];
    end
end
if isempty(eventDataFiles)
   eventObjects = [];
   return
end

% Identify unique event types
eventFileNames = {eventDataFiles.name};
eventTypes = regexp(eventFileNames, '(?<=event_data_).*(?=\.csv)', 'match', 'once');
eventTypeList = unique(eventTypes);

% Create event objects for each type and load data
disp('------------------------------------------');
disp('Loading event data files...')
eventObjects = struct();
for iType = 1:numel(eventTypeList)
    currType = eventTypeList{iType};
    currTypeFiles = eventFileNames(strcmp(eventTypes, currType));
    
    % Create an event object of the appropriate type
    switch currType
        case 'odor'
            eventObjects.(currType) = odorEvent();
        case 'optostim'
            eventObjects.(currType) = optoStimEvent();
        case 'panelsflash'
            eventObjects.(currType) = panelsFlashEvent();
        case 'soundstim'
            eventObjects.(currType) = soundStimEvent();
        case 'ir-laser'
            continue % I don't think I'm going to want to analyze these
        case {'locomotion', 'isolatedmovement', 'grooming', 'quiescence', 'flailing', 'ballstop'}
            eventObjects.(currType) = behaviorEvent(currType);
    end%switch
    
    % Import data from all event files of the current type
    for iFile = 1:numel(currTypeFiles)
        disp(['Loading ', currTypeFiles{iFile}, '...']);
        eventObjects.(currType) = eventObjects.(currType).load_csv(parentDir, 'fileNamePrefix', ...
            currTypeFiles{iFile}(1:10));
    end    
end
disp('All event data loaded')

end