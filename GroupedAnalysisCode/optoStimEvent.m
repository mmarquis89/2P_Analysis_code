classdef optoStimEvent < stimEvent

    
    properties (SetAccess = immutable)
        
    end
    
    methods
        
        % Constructor
        function obj = optoStimEvent(expID)
           obj = obj@stimEvent(expID, 'OptoStim');            
        end
    end
    
    
end