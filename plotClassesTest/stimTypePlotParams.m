classdef stimTypePlotParams < plotParams
        
    methods       
        % Constructor
        function obj = stimTypePlotParams(goodTrials, stimNames, stimTrialGroups)
           obj@plotParams(goodTrials);
           obj.groupType = 'StimTypes';
           obj.groupNames = stimNames;
           obj.trialGroups = goodTrials .* stimTrialGroups;       
        end
    end
    
    
    
end