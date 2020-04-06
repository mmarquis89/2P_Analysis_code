classdef odorEvent < stimEvent

    
    properties (SetAccess = immutable)
        odorName
        concentration
    end
    
    methods
        
        % Constructor
        function obj = odorEvent(expID, odorName, odorConc)
           obj = obj@stimEvent(expID, 'Odor');            
           obj.odorName = odorName;
           obj.concentration = odorConc;
        end
    end
    
    
end