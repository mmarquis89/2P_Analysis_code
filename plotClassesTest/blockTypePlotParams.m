classdef blockTypePlotParams < plotParams
    
    properties
        blockBounds
    end
    
    methods       
        % Constructor
        function obj = blockTypePlotParams(goodTrials, blockBounds, blockTypeNames, blockShading)
           obj@plotParams(goodTrials);
           obj.groupType = 'BlockTypes';
           obj.groupNames = blockTypeNames;
           obj.groupShading = blockShading;
           
           obj.blockBounds = blockBounds;
           
           % Divide trials into trial groups
           trialGroups = zeros(size(goodTrials));
           for iBound = 1:numel(blockBounds) - 1
               trialGroups(blockBounds(iBound)+1:blockBounds(iBound + 1)) = iBound;
           end
           trialGroups(groupBounds(end)+1:end) = iBound + 1;
           trialGroups(logical(mod(trialGroups, 2))) = 1;
           trialGroups(~logical(mod(trialGroups, 2))) = 2;
           obj.trialGroups = goodTrials .* trialGroups;
           obj.fileNameSuffix = '_Block_Type_Groups';
        end
    end
    
    
    
end