classdef odorEvent < stimEvent
% ==================================================================================================   
%    
% Properties:
%       metadataFieldNames
%
% Methods:
%       .append_shorthand(trialNums, stimTiming, trialDuration, odorName, odorConcentration,
%                         flowRate)
%
% ================================================================================================== 
    
    properties (SetAccess = immutable)
        metadataFieldNames
    end
    
    methods
        
        % Constructor
        function obj = odorEvent(expID)
           obj = obj@stimEvent(expID, 'Odor');            
           obj.metadataFieldNames = {'odorName', 'concentration', 'flowRate'};
        end
        
        % Modified append_shorthand function
        function obj = append_shorthand(obj, trialNums, stimTiming, trialDuration, odorName, ...
                    odorConcentration, flowRate)
            mdFieldNames = obj.metadataFieldNames;
            mdFieldValues = {odorName, odorConcentration, flowRate};
            obj = append_shorthand@stimEvent(obj, trialNums, stimTiming, trialDuration, ...
                    'MetadataFieldNames', mdFieldNames, 'MetadataFieldValues', mdFieldValues);
        end
        
    end%Methods    
end%Class 
