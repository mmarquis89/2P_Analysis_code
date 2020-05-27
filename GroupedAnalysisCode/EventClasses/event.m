classdef event
% ==================================================================================================   
%    
% Properties:
%       eventData
%       type                (immutable)
%       nExps               (dependent)
%       fieldNames          (dependent)
%       trialList           (dependent)
%
% Methods:
%       .append_data(expID, trialNums, eventTimes, metadataFieldNames, metadataFieldValues)
%       .export_csv(savePath, 'fileNamePrefix', 'fileNameSuffix')
%       .load_csv(parentDir, 'fileNamePrefix, 'fileNameSuffix')
%       .create_logical_array(nSamples, sampRate)
%
%       Also supports direct table indexing syntax, which is redirected to obj.eventData:
%           obj.(eventDataFieldName)
%           obj(1, 1:5) or obj(1, {'var1', 'var2'})
%           obj{...}
%
% Subclasses:
%       stimEvent
%           --> odorEvent, panelsFlashEvent, optoStimEvent, soundStimEvent
%       behaviorEvent
%
%
% ==================================================================================================

    properties
        eventData       % Table containing all the onset and offset times for each trial
        eventDataSubset % Copy of event data table subsetted by one or more metadata fields
    end
    properties (SetAccess = immutable)
        type             % The name of the event (i.e. 'odor', 'optoStim', 'locomotion', etc.)
    end
    properties (Dependent)
        % Calculated values
        nExps            % Number of distinct experiment IDs in the data   
        fieldNames       % Names of all fields in event data table, in view-friendly format
        trialList        % Table with all the unique expID - trialNum pairs in the eventData
    end
    
    methods
        
        % CONSTRUCTOR
        function obj = event(type)
            obj.type = type;
            obj.eventData = [];            
        end
        
        % GET/SET METHODS FOR DEPENDENT PROPERTIES
        function nTrials = get.nExps(obj)
            nTrials = numel(unique(obj.eventData.expID));
        end
        function fieldNames = get.fieldNames(obj)
            variableNames = obj.eventData.Properties.VariableNames';
            varNum = (1:numel(variableNames))';
            fieldNames = (table(varNum, variableNames));
        end
        function trialList = get.trialList(obj) 
            if ~isempty(obj.eventData)
                trialList = unique(obj.eventData(:, {'expID', 'trialNum'}));
            else
                trialList = [];
            end
        end
        
        % SUBSREF OVERLOAD FOR DIRECT TABLE INDEXING
        function sref = subsref(obj, s)
            if ~isempty(obj.eventData)
                switch s(1).type
                    case '.'
                        % Determine whether user is refering to obj properties or obj.eventData fields
                        if ismember(s(1).subs, obj.fieldNames.variableNames)
                            sref = builtin('subsref', obj.eventData, s);                            
                        else
                            sref = builtin('subsref', obj, s);
                        end
                    otherwise
                        % Other indexing types passed directly on to eventData table
                        sref = builtin('subsref', obj.eventData, s);
                end
            else
                % Don't even try to overload builtin if eventData is empty
                sref = builtin('subsref', obj, s);
            end
        end
        
        % APPEND DATA FOR NEW TRIALS
        function obj = append_data(obj, expID, trialNums, eventTimes, metadataFieldNames, ...
                    metadataFieldValues)
            
            % Make sure input arguments are the right size
            if numel(trialNums) ~= numel(eventTimes)
                error('ERROR: eventTimes must contain one cell for each trial in trialNums')
            end
            
            % Create new table to append to eventData
            tbAppend = [];
            for iTrial = 1:numel(trialNums)
                newRow = table({expID}, trialNums(iTrial), {eventTimes{iTrial}(:, 1)}, ...
                            {eventTimes{iTrial}(:, 2)}, 'VariableNames', {'expID', 'trialNum', ...
                            'onsetTimes', 'offsetTimes'});
                tbAppend = [tbAppend; newRow];                
            end
            tbAppendFlat = flatten_table(tbAppend);
            
            % Add any extra metadata fields if necessary
            if ~isempty(metadataFieldNames)
                for iField = 1:numel(metadataFieldNames)
                    currMdVal = metadataFieldValues{iField};
                    if ~isscalar(currMdVal) || ~isnumeric(currMdVal)
                        currMdVal = {currMdVal};
                    end
                    tbAppendFlat.(metadataFieldNames{iField}) = repmat(currMdVal, ...
                        size(tbAppendFlat, 1), 1);
                end
            end
            
            % Flatten data into [trialNum, onsetTime, offsetTime] format
            obj.eventData = [obj.eventData; tbAppendFlat];
            
            % Delete any duplicate rows that were introduced
            obj.eventData = unique(obj.eventData);
        end
        
        % SAVE EVENT TIMES TO .CSV FILE
        function obj = export_csv(obj, savePath, varargin)
            
            [fileNamePrefix, fileNameSuffix] = process_filename_args(varargin);
            
            saveFileName = [fileNamePrefix 'event_data_', lower(obj.type), fileNameSuffix, '.csv'];
            
            % Write to file
            overwrite = 1;
            if exist(fullfile(savePath, saveFileName), 'file') ~=0
                dlgAns = questdlg(['Saving this file will overwrite an existing file in this', ...
                        ' directory...are you sure you want to proceed?'], ...
                        'Warning', 'Yes', 'No', 'No');
                if strcmp(dlgAns, 'No')
                   overwrite = 0;
                   disp('CSV export cancelled')
                end
            end
            if overwrite
                writetable(obj.eventData, fullfile(savePath, saveFileName));
            end
        end
         
        % IMPORT DATA FROM .CSV FILE
        function obj = load_csv(obj, parentDir, varargin)
            
            [fileNamePrefix, fileNameSuffix] = process_filename_args(varargin);        
            
            % Load file 
            tbAppend = readtable(fullfile(parentDir, [fileNamePrefix, 'event_data_', ...
                    lower(obj.type), fileNameSuffix, '.csv']));
                                       
            % Append to eventData table
            obj.eventData = [obj.eventData; tbAppend];
            
            % Delete any duplicate rows that were introduced
            if any(size(obj.eventData) ~= size(unique(obj.eventData)))
                obj.eventData = unique(obj.eventData);
                warning('duplicate rows removed from event data after import');
            end
            
        end
        
        % CONVERT EVENT TIME INTO LOGICAL ARRAY 
        function outputArr = create_logical_array(obj, nSamples, sampTimes, trialList)
            % TrialList = table with columns [expID][trialNum][onsetTime][offsetTime]
            if nargin < 4
                trialList = obj.trialList;
            end
            if numel(sampTimes) ~= nSamples
               error('Number of elements in sampTimes must equal nSamples') 
            end
            nTrials = size(trialList, 1);
                
            outputArr = zeros(nSamples, nTrials); % --> [sample, trial]
            for iTrial = 1:nTrials
                currTrialEvents = innerjoin(trialList(iTrial, :), obj.eventData);
                onsetSamples = [];
                offsetSamples = [];
                for iEvent = 1:size(currTrialEvents, 1)
                    [~, onsetSamples(iEvent)] = min(abs(volTimes - ...
                            currTrialEvents.onsetTime(iEvent)));
                    [~, offsetSamples(iEvent)] = min(abs(volTimes - ...
                            currTrialEvents.offsetTime(iEvent)));
                end
                onsetSamples(onsetSamples < 1) = 1;
                offsetSamples(offsetSamples > nSamples) = nSamples;
                for iEvent = 1:numel(onsetSamples)
                    outputArr(onsetSamples(iEvent):offsetSamples(iEvent), iTrial) = 1;
                end
            end            
        end
        
    end%Methods
    
end%Class

% ==================================================================================================
% Local functions
% ==================================================================================================

% Flatten table
function tbOut = flatten_table(tbIn)
tbOut = [];
for iTrial = 1:size(tbIn, 1)
    currExpID = tbIn.expID{iTrial};
    currTrialNum = tbIn.trialNum(iTrial);
    currOnsets = tbIn.onsetTimes{iTrial};
    currOffsets = tbIn.offsetTimes{iTrial};
    for iEvent = 1:numel(currOnsets)
        newRow = table({currExpID}, currTrialNum, currOnsets(iEvent), currOffsets(iEvent), ...
                'VariableNames', {'expID', 'trialNum', 'onsetTime', 'offsetTime'});
        tbOut = [tbOut; newRow];
    end
end
end 

% Process optional arguments for .csv read/write functions
function [fileNamePrefix, fileNameSuffix] = process_filename_args(varargin)

p = inputParser;
addParameter(p, 'fileNamePrefix', '');
addParameter(p, 'fileNameSuffix', '');
parse(p, varargin{1}{:});
fileNamePrefix = p.Results.fileNamePrefix;
fileNameSuffix = p.Results.fileNameSuffix;

% Make sure suffix (if provided) starts with underscore
if ~isempty(fileNameSuffix) && ~strcmp(fileNameSuffix(1), '_')
    fileNameSuffix = ['_', fileNameSuffix];
end

% Make sure prefix (if provided) ends with underscore
if ~isempty(fileNamePrefix) && ~strcmp(fileNamePrefix(end), '_')
    fileNamePrefix = [fileNamePrefix, '_'];
end

end




