classdef plotParams
    
   properties
       
      % Group name
      groupType
      
      % Will vary across subclass instances, but all will have these properties
      trialGroups
      groupNames
      groupShading
      fileNameSuffix
      plotTitleSuffix
      
      % Common across all plots
      shadeEpochs
      shadeEpochColors
      shadeEpochNames
      
   end
   
   methods
       % Constructor function
       function obj = plotParams(goodTrials)
           obj.trialGroups = goodTrials;
       end
       
       % Add stim shading epochs
       function obj = add_shadeEpochs(obj, epochTimes, epochColors, epochNames)
           obj.shadeEpochs = epochTimes;
           obj.shadeEpochColors = epochColors;
           obj.shadeEpochNames = epochNames;
       end
       
       % Generate plot title suffix
       function obj = get_plotSuffix(obj)
          obj.plotTitleSuffix = make_plotTitleSuffix(obj.groupNames); 
       end
       
       % Generate file name suffix
       function obj = get_fileNameSuffix(obj)
           obj.fileNameSuffix = make_fileNameSuffix(obj.groupNames); 
       end
       
   end  
    
end