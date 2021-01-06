function split_roiDefs(parentDir, splitRois)
%===================================================================================================
% 
% Splits an ROI and makes of its constituent subROIs a new ROI with a number appended to the 
% name, and re-saves it to the roiDefs. Keeps the original ROI with no number appended to the name. 
% 
% INPUTS: 
%       parentDir   = a directory containing the roiDef files you want to process (file name must 
%                     start with 'roiDefs')
%       splitRois   = a cell array containing the names of the ROIs you want to split (pass an 
%                     empty vector to split all ROIs)
%
%===================================================================================================

roiDefFiles = dir(fullfile(parentDir, ['roiDefs_trial*.mat']));

for iFile = 1:numel(roiDefFiles)
    
    % Load roiDef data
    load(fullfile(parentDir, roiDefFiles(iFile).name), 'roiDefs');
    
    % Split Rois 
    if isempty(splitRois)
        splitRois = {roiDefs.name};
    end
    newRoiDefs = roiDefs;
    for iRoi = 1:numel(splitRois)
        currRoiName = splitRois{iRoi};
        subRois = roiDefs(strcmp({roiDefs.name}, currRoiName)).subROIs;
        if numel(subRois) > 1
            for iSubRoi = 1:numel(subRois)
                newRoiDef = struct();
                newRoiDef.name = [currRoiName, '-', num2str(iSubRoi)];
                newRoiDef.subROIs = subRois(iSubRoi);
                newRoiDef.color = roiDefs(strcmp({roiDefs.name}, currRoiName)).color;
                newRoiDefs = [newRoiDefs, newRoiDef];
            end
        end
    end
    
    % Save modified roiDefs file
    roiDefs = newRoiDefs;
    save(fullfile(parentDir, roiDefFiles(iFile).name), 'roiDefs');
        
end

end