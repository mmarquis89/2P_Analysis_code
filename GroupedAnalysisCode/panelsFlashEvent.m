classdef panelsFlashEvent < stimEvent

    
    properties (SetAccess = immutable)
        
    end
    
    methods
        
        % Constructor
        function obj = visualStimEvent(expID, brightness)
           obj = obj@stimEvent(expID, 'PanelsFlash');            
        end
    end
    
    
end