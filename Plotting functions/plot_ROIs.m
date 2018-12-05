function plot_ROIs(ROImetadata)

% Create a figure for each ROI and plot the ROI in each plane on a separate subplot
for iROI = 1:numel(ROImetadata)
   figure(iROI); clf
   nSubplots = numel(ROImetadata{iROI});
   nPlots = numSubplots(nSubplots);
   for iSubplot = 1:nSubplots
       
      % Plot reference image 
      subaxis(nPlots(1), nPlots(2), iSubplot); 
      imshow(ROImetadata{iROI}(iSubplot).refImg, [])
%       imagesc(ROImetadata{iROI}(iSubplot).refImg);
      % Overlay ROI
      hold on
      currData = ROImetadata{iROI}(iSubplot);
      plot(currData.xi, currData.yi, 'LineWidth', 2)
      title(['Plane #', num2str(currData.plane)])
      
   end
   
   % Add title for entire figure
   suptitle(['ROI #', num2str(iROI)])
   
end%for

end%function