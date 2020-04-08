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
            mdFieldNames = p.Results.MetadataFieldNames;
            mdFieldValues = p.Results.MetadataFieldValues;
            
            % Make sure any extra field name/value pairs match up 
            if numel(mdFieldNames) ~= numel(mdFieldValues)
                error('ERROR: MetadataFieldNames and MetadataFieldValues must be same length');
            end
            
            % Generate stim times from shorthand input
            eventTimes = {};
            for iTrial = 1:numel(trialNums)
               onsetTimes = stimTiming(1):sum(stimTiming(2:3)):trialDuration;
               offsetTimes = onsetTimes + stimTiming(2);
               eventTimes{iTrial} = [onsetTimes', offsetTimes'];
            end           
            
            obj = obj.append_data(trialNums, eventTimes, mdFieldNames, mdFieldValues);
        end
    end
    
    
end