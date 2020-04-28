classdef soundStimEvent < stimEvent
% ==================================================================================================   
%    
% Properties:
%       metadataFieldNames
%
% Methods:
%       .append_shorthand(expID, trialNums, stimTiming, trialDuration, toneFreq, ampSetting)
%
% ================================================================================================== 
    properties (SetAccess = immutable)
        metadataFieldNames
    end
    
    methods
        
        % Constructor
        function obj = soundStimEvent()
            obj = obj@stimEvent('SoundStim');
            obj.metadataFieldNames = {'toneFreq', 'ampSetting'};
        end
        
        % Modified append_shorthand function
        function obj = append_shorthand(obj, expID, trialNums, stimTiming, trialDuration, toneFreq, ...
                    ampSetting)
            mdFieldNames = obj.metadataFieldNames;
            mdFieldValues = {toneFreq, ampSetting};
            obj = append_shorthand@stimEvent(obj, expID, trialNums, stimTiming, trialDuration, ...
                    'MetadataFieldNames', mdFieldNames, 'MetadataFieldValues', mdFieldValues);
        end
        
    end
    
    
end