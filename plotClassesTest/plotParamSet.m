classdef plotParams
    
   properties
      % Common across all plots
      shadeEpochs
      shadeEpochColors
      shadeEpochNames
       
      % Will vary across subclass instances, but all will have these properties
      trialGroups
      groupNames
      fileNameSuffix
      plotTitleSuffix
      
   end
   
   methods
      
       function obj = plotParams(obj, groupNames)
           
       end
       
   end  
    
end