classdef stimEvent < event

    
    properties
        
    end
    
    methods
        
        % Constructor
        function obj = stimEvent(expID, type)
           obj = obj@event(expID, type);            
        end
        
        % Append using [pre-stim, duration, ISI] notation
        function obj = append_shorthand(obj, trialNums, stimTiming, trialDuration, varargin)
            
            % Parse optional arguments
            p = inputParser;
            addParameter(p, 'MetadataFieldNames', []);
            addParameter(p, 'MetadataFieldValues', []);
            parse(p, varargin{:});
            extraFieldNames = p.Results.MetadataFieldNames;
            extraFieldValues = p.Results.MetadataFieldValues;
            
            % Make sure any extra field name/value pairs match up 
            if numel(extraFieldNames) ~= numel(extraFieldValues)
                error('ERROR: MetadataFieldNames and MetadataFieldValues must be same length');
            end
            
            % Generate stim times from shorthand input
            eventTimes = {};
            for iTrial = 1:numel(trialNums)
               onsetTimes = stimTiming(1):sum(stimTiming(2:3)):trialDuration;
               offsetTimes = onsetTimes + stimTiming(2);
               eventTimes{iTrial} = [onsetTimes', offsetTimes'];
            end           
            
            obj = obj.append_data(trialNums, eventTimes, extraFieldNames, extraFieldValues);
        end
    end
    
    
end