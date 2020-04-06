classdef behaviorEvent < event

    
    properties
        
    end
    
    methods
        
        % Constructor
        function obj = behaviorEvent(expID, type)
           obj = obj@event(expID, type);            
        end
        
    end
    
    
end