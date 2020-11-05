function combine_roiDefs(parentDir, combROIs, newROI)
%===================================================================================================
% 
% Combines two or more ROI defs into a single, new ROI for all roiDef files in a given directory. 
% Then adds the new ROI to the roiDefs file and re-saves it. Will skip any files that already 
% contain an ROI with the same name as newROI.
% 
% INPUTS: 
%       parentDir   = a directory containing the roiDef files you want to process (file name must 
%                     start with 'roiDefs')
%       combROIs    = a cell array containing the names of the ROIs you want to group together
%       newROI      = the name for the newly created ROI
%
%===================================================================================================

roiDefFiles = dir(fullfile(parentDir, ['roiDefs*.mat']));

for iFile = 1:numel(roiDefFiles)
    
    % Load roiDef data
    load(fullfile(parentDir, roiDefFiles(iFile).name), 'roiDefs');
    
    % Create combined ROI 
    newRoiDef = struct();
    newRoiDef.name = newROI;
    newRoiDef.subROIs = [roiDefs(ismember({roiDefs.name}, combROIs)).subROIs];
    if ~isempty(newRoiDef.subROIs) && ~any(strcmp({roiDefs.name}, newROI))
        newRoiDef.color = roiDefs(1).color;
        
        % Add to original ROI defs
        roiDefs = [roiDefs, newRoiDef];
        
        % Save modified roiDefs file
        save(fullfile(parentDir, roiDefFiles(iFile).name), 'roiDefs');
    end
end

end