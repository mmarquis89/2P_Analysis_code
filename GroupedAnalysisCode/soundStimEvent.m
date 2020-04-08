classdef soundStimEvent < stimEvent

    properties (SetAccess = immutable)
        metadataFieldNames
    end
    
    methods
        
        % Constructor
        function obj = soundStimEvent(expID, ampSetting)
           obj = obj@stimEvent(expID, 'SoundStim');            
           obj.ampSetting = ampSetting;
        end
        
    end
    
    
end