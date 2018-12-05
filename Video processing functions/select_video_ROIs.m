function roiData = select_video_ROIs(parentDir)
%===================================================================================================
% DEFINE ROI FOR FLY MOVEMENT IN BEHAVIOR VIDEO 
%
% Prompts the user to draw an ROI on a frame of the behavior video for the purposes of analyzing 
% optic flow. The ROI should be drawn around the area where fly movement is occuring, ideally 
% avoiding areas that the ball might enter if it shakes/bounces. Saves the ROI mask for future use 
% in creation of combined optic flow + behavior videos.
%
% INPUTS:
%   parentDir = the directory containing the snapshot of the behavior video
%
%===================================================================================================

if nargin < 1
   parentDir = ''; 
end

[imgFile, imgDir, ~] = uigetfile(fullfile(parentDir, '*'));

% Read image (or extract first frame if image is a video)
if regexp(imgFile, '.*(.avi|.mp4|.mov)')
    myVid = VideoReader(fullfile(imgDir, imgFile));
    firstFrame = readFrame(myVid);
    plotImg = uint8(firstFrame(:,:,1));
else
    plotImg = uint8(imread(fullfile(imgDir, imgFile)));
end

h = figure(1);clf; hold on
im = imshow(plotImg, []);

% Prompt user to create an ROI
h.Name = 'Define an ROI for fly movement';
disp('Define an ROI for fly movement')
roiData = []; xi = []; yi = [];
[roiData, ~, ~] = roipoly;

% Save ROI data
saveDir = uigetdir(parentDir, 'Select a save directory');
if saveDir == 0
    % Throw error if user canceled without choosing a directory
    disp('ERROR: you must select a save directory or provide one as an argument');
else
    % Prompt user for file name
    fileName = inputdlg('Please choose a file name', 'Save ROI data', 1, {'Behavior_Vid_ROI_Data'});
    fileName = fileName{:};
    
    % Save data
    save(fullfile(saveDir, [fileName, '.mat']), 'roiData', '-v7.3');
end
disp('ROI data saved')

end