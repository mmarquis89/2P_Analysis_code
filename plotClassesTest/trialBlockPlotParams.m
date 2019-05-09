classdef trialBlockPlotParams < plotParams
    
    properties
        blockBounds
    end
    
    methods       
        % Constructor
        function obj = trialBlockPlotParams(goodTrials, blockBounds, blockNames, blockShading)
           obj@plotParams(goodTrials);
           obj.groupType = 'TrialBlocks';
           obj.groupNames = blockNames;
           obj.groupShading = blockShading;
           
           obj.blockBounds = blockBounds;
           
           % Define trial groups
           trialGroups = zeros(size(goodTrials));
           for iBound = 1:numel(blockBounds) - 1
               trialGroups(blockBounds(iBound)+1:blockBounds(iBound + 1)) = iBound;
           end
           trialGroups(groupBounds(end)+1:end) = iBound + 1;
           obj.trialGroups = goodTrials .* trialGroups;
           obj.fileNameSuffix = '_Blocks_Separated';
        end
    end
    
    
    
end