classdef panelsFlashEvent < stimEvent
    
    properties (SetAccess = immutable)
        metadataFieldNames
    end
    
    methods
        
        % Constructor
        function obj = panelsFlashEvent()
           obj = obj@stimEvent('PanelsFlash'); 
           obj.metadataFieldNames = {'brightness', 'type'};
        end
        
        % Modified append_data function
        function obj = append_data(obj, expID, trialNums, eventTimes, brightness, patternType)
            mdFieldNames = obj.metadataFieldNames;
            mdFieldValues = {brightness, patternType};
            obj = append_data@stimEvent(obj, expID, trialNums, eventTimes, ...
                    mdFieldNames, mdFieldValues);
        end
    end
end