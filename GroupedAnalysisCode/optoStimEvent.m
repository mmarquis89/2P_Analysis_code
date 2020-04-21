classdef optoStimEvent < stimEvent
% ==================================================================================================   
%    
% Properties:
%       metadataFieldNames
%
% Methods:
%       .append_shorthand(trialNums, stimTiming, trialDuration, LEDpower, dutyCycle)
%
% ================================================================================================== 
    
    properties (SetAccess = immutable)
        metadataFieldNames
    end
    
    methods
        
        % Constructor
        function obj = optoStimEvent()
            obj = obj@stimEvent('OptoStim');
            obj.metadataFieldNames = {'LEDpower', 'dutyCycle'};
        end
        
        % Modified append_shorthand function
        function obj = append_shorthand(obj, expID, trialNums, stimTiming, trialDuration, LEDpower, ...
                dutyCycle)
            mdFieldNames = obj.metadataFieldNames;
            mdFieldValues = {LEDpower, dutyCycle};
            obj = append_shorthand@stimEvent(obj, expID, trialNums, stimTiming, trialDuration, ...
                'MetadataFieldNames', mdFieldNames, 'MetadataFieldValues', mdFieldValues);
        end
    end
    
    
end