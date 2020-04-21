classdef event
% ==================================================================================================   
%    
% Properties:
%       eventData
%       eventDataSubset
%       type                (immutable)
%       nExps               (dependent)
%       fieldNames          (dependent)
%
% Methods:
%       .append_data(expID, trialNums, eventTimes, metadataFieldNames, metadataFieldValues)
%       .export_csv(savePath, 'fileNamePrefix', 'fileNameSuffix')
%       .import_csv(parentDir, 'fileNamePrefix, 'fileNameSuffix')
%       .metadata_subset(fieldName, values)
%       .create_logical_array(nSamples, sampRate)
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
                newRow = table(expID, trialNums(iTrial), {eventTimes{iTrial}(:, 1)}, ...
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
        function export_csv(obj, savePath, varargin)
            
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
        
        % SUBSET EVENTS BY TRIAL NUMBER
        function obj = trial_subset(obj, trialNums)
            obj.eventData = obj.eventData(ismember(obj.eventData.trialNum, trialNums), :);
        end
        
        % SUBSET EVENTS BY EXTRA METADATA FIELD VALUES
        function obj = metadata_subset(obj, fieldName, values)
            obj.eventDataSubset = obj.eventData(ismember(obj.eventData.(fieldName), values), :);
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
        newRow = table(currExpID, currTrialNum, currOnsets(iEvent), currOffsets(iEvent), ...
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
parse(p, varargin{:});
fileNamePrefix = p.fileNamePrefix;
fileNameSuffix = p.fileNameSuffix;

% Make sure suffix (if provided) starts with underscore
if ~isempty(fileNameSuffix) && ~strcmp(fileNameSuffix(1), '_')
    fileNameSuffix = ['_', fileNameSuffix];
end

% Make sure prefix (if provided) ends with underscore
if ~isempty(fileNamePrefix) && ~strcmp(fileNamePrefix(end), '_')
    fileNamePrefix = [fileNamePrefix, '_'];
end

end




