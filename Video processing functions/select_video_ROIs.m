function select_video_ROIs(parentDir, sid)
%===================================================================================================
% DEFINE ROI FOR FLY MOVEMENT IN BEHAVIOR VIDEO 
%
% Prompts the user to draw an ROI on a frame of the behavior video for the purposes of analyzing 
% optic flow. The ROI should be drawn around the area where fly movement is occuring, ideally 
% avoiding areas that the ball might enter if it shakes/bounces. Saves the ROI mask for future use 
% in creation of combined optic flow + behavior videos.
%
% INPUTS:
%   parentDir = the directory containing the behavior videos for the experiment you want to analyze
%
%   sid       = the session ID that you want to define an ROI for
%===================================================================================================

% Load and plot 1st frame acquired during the first trial (number choice is arbitrary)
myFolders = dir(fullfile(parentDir, ['*sid_', num2str(sid), '*tid*']));
myFolders = myFolders([myFolders.isdir]);
if isempty(myFolders)
    myFolders = dir(fullfile(parentDir, ['*sid_', num2str(sid), '*bid*']));
    myFolders = myFolders([myFolders.isdir]);
end

iTrial = 1;
myFrames = dir(fullfile(parentDir, myFolders(iTrial).name, '*.tif'));
while isempty(myFrames) % In case the first trial doesn't have any video frames
    iTrial = iTrial + 1;
    myFrames = dir(fullfile(parentDir, myFolders(iTrial).name, '*.tif'));
end
plotImg = uint8(imread(fullfile(parentDir, myFolders(iTrial).name, myFrames(1).name)));
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