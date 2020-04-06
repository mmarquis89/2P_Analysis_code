classdef event
% ==================================================================================================   
%    
% Properties:
%       expID               (immutable)
%       eventTimes
%       type                (immutable)
%       nTrials             (dependent)
%
% Methods:
%       .append_data(trialNums, eventTimes)
%       .export_csv(savePath, fileNameSuffix)
%       .import_csv(parentDir, fileNameSuffix)
%       .subset_trials(trialNums)
%
%
% Subclasses:
%       odorEvent
%
%
%
% ==================================================================================================

    properties
        eventTimes       % Table containing all the onset and offset times for each trial
    end
    properties (SetAccess = immutable)
        expID            % ID of the experiment the events were ocurring in
        type             % The name of the event (i.e. 'odor', 'optoStim', 'locomotion', etc.)
    end
    properties (Dependent)
        % Calculated values
        nTrials          % Number of distinct trials (with acquisition gaps) in the data   
    end
    
    methods
        
        % CONSTRUCTOR
        function obj = event(expID, type)
            obj.expID = expID;
            obj.type = type;
            obj.eventTimes = [];            
        end
        
        % GET/SET METHODS FOR DEPENDENT PROPERTIES
        function nTrials = get.nTrials(obj)
            nTrials = numel(unique(obj.eventTimes.trialNum));
        end
        
        % APPEND DATA FOR NEW TRIALS
        function obj = append_data(obj, trialNums, eventTimes)
            
            % Make sure input arguments are the right size
            if numel(trialNums) ~= numel(eventTimes)
                error('ERROR: eventTimes must contain one cell for each trial in trialNums')
            end
            
            % Make sure data doesn't already exist for those trial numbers
            if ~isempty(obj.eventTimes) && ...
                        ~isempty(intersect(obj.eventTimes.trialNum, trialNums))
               error('ERROR: at least one of those trial numbers already exists in the event data')
            end
            
            % Append new data to eventTimes table
            tbAppend = [];
            for iTrial = 1:numel(trialNums)
                newRow = table(trialNums(iTrial), {eventTimes{iTrial}(:, 1)}, ...
                            {eventTimes{iTrial}(:, 2)}, 'VariableNames', {'trialNum', ...
                            'onsetTimes', 'offsetTimes'});
                tbAppend = [tbAppend; newRow];                
            end
            
            % Flatten data into [trialNum, onsetTime, offsetTime] format
            obj.eventTimes = [obj.eventTimes; flatten_table(tbAppend)];
        end
        
        % SAVE EVENT TIMES TO .CSV FILE
        function export_csv(obj, savePath, fileNameSuffix)
            
            % Make sure suffix starts with underscore
            if ~isempty(fileNameSuffix) && ~strcmp(fileNameSuffix(1), '_') 
               fileNameSuffix = ['_', fileNameSuffix]; 
            end
            
            % Write to file
            writetable(obj.eventTimes, fullfile(savePath, ['event_data_', lower(obj.type), '_', ...
                    obj.expID, fileNameSuffix, '.csv']));
        end
         
        % IMPORT DATA FROM .CSV FILE
        function obj = load_csv(obj, parentDir, fileNameSuffix)
            
            % Make sure suffix starts with underscore
            if ~isempty(fileNameSuffix) && ~strcmp(fileNameSuffix(1), '_')
                fileNameSuffix = ['_', fileNameSuffix];
            end
            
            % Load file 
            tbAppend = csvread(fullfile(parentDir, ['event_data_', lower(obj.type), '_', obj.expID, ...
                    fileNameSuffix, '.csv']), 1, 0);
                            
            % Make sure data doesn't already exist for those trial numbers
            if ~isempty(obj.eventTimes) && ...
                        ~isempty(intersect(obj.eventTimes.trialNum, tbAppend.trialNum))
               error('ERROR: at least one of those trial numbers already exists in the event data')
            end
            
            % Append to eventTimes table
            obj.eventTimes = [obj.eventTimes; tbAppend];
        end
        
        % SELECT SUBSET OF ALL TRIALS
        function obj = subset_trials(obj, trialNums)
            obj.eventTimes = obj.eventTimes(ismember(obj.eventTimes.trialNum, trialNums), :);
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

