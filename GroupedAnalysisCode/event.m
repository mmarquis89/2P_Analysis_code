classdef event
% ==================================================================================================   
%    
% Properties:
%       expID               (immutable)
%       eventData
%       type                (immutable)
%       nTrials             (dependent)
%       fieldNames          (dependent)
%
% Methods:
%       .append_data(trialNums, eventTimes, metadataFieldNames, metadataFieldValues)
%       .clear_data()
%       .export_csv(savePath, fileNameSuffix)
%       .import_csv(parentDir, fileNameSuffix)
%       .trial_subset(trialNums)
%       .metadata_subset(fieldName, values)
%
% Subclasses:
%       stimEvent
%       flailingEvent
%
%
% ==================================================================================================

    properties
        eventData       % Table containing all the onset and offset times for each trial
    end
    properties (SetAccess = immutable)
        expID            % ID of the experiment the events were ocurring in
        type             % The name of the event (i.e. 'odor', 'optoStim', 'locomotion', etc.)
    end
    properties (Dependent)
        % Calculated values
        nTrials          % Number of distinct trials (with acquisition gaps) in the data   
        fieldNames       % Names of all fields in event data table, in view-friendly format
    end
    
    methods
        
        % CONSTRUCTOR
        function obj = event(expID, type)
            obj.expID = expID;
            obj.type = type;
            obj.eventData = [];            
        end
        
        % GET/SET METHODS FOR DEPENDENT PROPERTIES
        function nTrials = get.nTrials(obj)
            nTrials = numel(unique(obj.eventData.trialNum));
        end
        function fieldNames = get.fieldNames(obj)
            variableNames = obj.eventData.Properties.VariableNames';
            varNum = (1:numel(variableNames))';
            fieldNames = (table(varNum, variableNames));
        end
        
        % APPEND DATA FOR NEW TRIALS
        function obj = append_data(obj, trialNums, eventTimes, metadataFieldNames, ...
                    metadataFieldValues)
            
            % Make sure input arguments are the right size
            if numel(trialNums) ~= numel(eventTimes)
                error('ERROR: eventTimes must contain one cell for each trial in trialNums')
            end
            
            % Warn user if data already exists for those trial numbers
            if ~isempty(obj.eventData) && ...
                        ~isempty(intersect(obj.eventData.trialNum, trialNums))
                warning(['at least one of those trial numbers already exists in the ', ...
                        'event data']);
            end
            
            % Create new table to append to eventData
            tbAppend = [];
            for iTrial = 1:numel(trialNums)
                newRow = table(trialNums(iTrial), {eventTimes{iTrial}(:, 1)}, ...
                            {eventTimes{iTrial}(:, 2)}, 'VariableNames', {'trialNum', ...
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
        
        % CLEAR ALL EVENT DATA
        function obj = clear_data(obj)
            obj.eventData = [];
        end 
        
        % SAVE EVENT TIMES TO .CSV FILE
        function export_csv(obj, savePath, fileNameSuffix)
            
            % Make sure suffix starts with underscore
            if ~isempty(fileNameSuffix) && ~strcmp(fileNameSuffix(1), '_') 
               fileNameSuffix = ['_', fileNameSuffix]; 
            end
            
            % Write to file
            writetable(obj.eventData, fullfile(savePath, [obj.expID, '_event_data_', ...
                    lower(obj.type), fileNameSuffix, '.csv']));
        end
         
        % IMPORT DATA FROM .CSV FILE
        function obj = load_csv(obj, parentDir, fileNameSuffix)
            
            % Make sure suffix starts with underscore
            if ~isempty(fileNameSuffix) && ~strcmp(fileNameSuffix(1), '_')
                fileNameSuffix = ['_', fileNameSuffix];
            end
            
            % Load file 
            tbAppend = csvread(fullfile(parentDir, [obj.expID, '_event_data_', lower(obj.type),  ...
                    fileNameSuffix, '.csv']), 1, 0);
                            
            % Make sure data doesn't already exist for those trial numbers
            if ~isempty(obj.eventData) && ...
                        ~isempty(intersect(obj.eventData.trialNum, tbAppend.trialNum))
                warning(['at least one of those trial numbers already exists in the ', ...
                        'event data']);
            end
            
            % Append to eventData table
            obj.eventData = [obj.eventData; tbAppend];
            
            % Delete any duplicate rows that were introduced
            obj.eventData = unique(obj.eventData);
        end
        
        % SUBSET EVENTS BY TRIAL NUMBER
        function obj = trial_subset(obj, trialNums)
            obj.eventData = obj.eventData(ismember(obj.eventData.trialNum, trialNums), :);
        end
        
        % SUBSET EVENTS BY EXTRA METADATA FIELD VALUES
        function obj = metadata_subset(obj, fieldName, values)
            obj.eventData = obj.eventData(ismember(obj.eventData.(fieldName), values), :);
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
                currTrialNum = tbIn.trialNum(iTrial);
                currOnsets = tbIn.onsetTimes{iTrial};
                currOffsets = tbIn.offsetTimes{iTrial};
                for iEvent = 1:numel(currOnsets)
                    newRow = table(currTrialNum, currOnsets(iEvent), currOffsets(iEvent), ...
                            'VariableNames', {'trialNum', 'onsetTime', 'offsetTime'});
                    tbOut = [tbOut; newRow];
                end
            end
end 

