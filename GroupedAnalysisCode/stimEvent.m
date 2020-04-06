classdef stimEvent < event

    
    properties
        
    end
    
    methods
        
        % Constructor
        function obj = stimEvent(expID, type)
           obj = obj@event(expID, type);            
        end
        
        % Append using [pre-stim, duration, ISI] notation
        function obj = append_shorthand(obj, trialNums, stimTiming, trialDuration)
            
            eventTimes = {};
            for iTrial = 1:numel(trialNums)
               onsetTimes = stimTiming(1):sum(stimTiming(2:3)):trialDuration;
               offsetTimes = onsetTimes + stimTiming(2);
               eventTimes{iTrial} = [onsetTimes', offsetTimes'];
            end
            
            obj = obj.append_data(trialNums, eventTimes);
        end
    end
    
    
end