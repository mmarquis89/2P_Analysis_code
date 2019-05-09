classdef alternatingTrialsPlotParams < plotParams
    
    properties
        baselineBound
        omitBaseline = 0;
    end
    
    methods       
        % Constructor
        function obj = alternatingTrialsPlotParams(goodTrials, ...
                    baselineBound, groupNames, groupShading)
           obj@plotParams(goodTrials);
           obj.groupType = 'AlternatingTrials';
           obj.groupNames = groupNames;
           obj.groupShading = groupShading;
           
           obj.blockBounds = blockBounds;
           
           % Divide trials into trial groups
           trialGroups = zeros(size(goodTrials));
           trialGroups(1:baselineBound-1) = 1; % Baseline period
           trialGroups(baselineBound:2:end) = 2;
           trialGroups(baselineBound+1:2:end) = 3;
           obj.trialGroups = goodTrials .* trialGroups;
           obj.fileNameSuffix = '_Alternating_Trials';
        end
        
        % Drop trials from the baseline block
        function obj = drop_baseline(obj)
           obj.groupNames = obj.groupNames(2:3); 
           obj.groupShading = obj.groupShading(2:3);
           obj.trialGroups = obj.trialGroups - 1;
           obj.trialGroups(obj.trialGroups < 0) = 0;
           obj.omitBaseline = 1;
           obj.fileNameSuffix = [obj.fileNameSuffix, '_No_Baseline'];
        end
        
        
    end
    
    
    
end