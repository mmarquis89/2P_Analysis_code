classdef panelsFlashEvent < stimEvent
    
    properties (SetAccess = immutable)
        metadataFieldNames
    end
    
    methods
        
        % Constructor
        function obj = panelsFlashEvent(expID)
           obj = obj@stimEvent(expID, 'PanelsFlash'); 
           obj.metadataFieldNames = {'brightness', 'type'};
        end
        
        % Modified append_data function
        function obj = append_data(obj, trialNums, eventTimes, brightness, patternType)
            mdFieldNames = obj.metadataFieldNames;
            mdFieldValues = {brightness, patternType};
            obj = append_data@stimEvent(obj, trialNums, eventTimes, ...
                    mdFieldNames, mdFieldValues);
        end
    end
end