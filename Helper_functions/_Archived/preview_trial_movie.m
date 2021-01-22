function preview_trial_movie(sessionObj, planeNum, trialNum, cMapRange, sigma, fps)
%===================================================================================================
% Opens the movie player app and loads a movie of a specific trial and plane (useful for screening
% trials for movement when deciding which trial and volume to use as a registration template). It 
% also enlarges the movie for easier viewing and smooths the data slightly (this can be deactivated
% by passing a very small (but nonzero) value for sigma).
%
% INPUTS:
%
%   sessionObj = a matfile object containing imaging data as a variable 'wholeSession' with 
%                dimensions [y, x, plane, volume, trial]
%
%   planeNum    = the number of the plane to show video from.
%
%   trialNum    = the number of the trial to show video from.
%
%   cMapRange   = <OPTIONAL> a 1x2 vector with the [min, max] values for the colormap scaling. Can 
%                pass a scalar value instead to use as the max map value and set the min value to 
%                zero. Pass [] to use the default range of [0 5000].
%
%   sigma       = <OPTIONAL> the standard deviation for the gaussian smoothing function. Pass [] to 
%                 use the default value of 0.5.
%
%   fps         = <OPTIONAL> the frame rate to play the video at. Pass [] to default to 30 fps.
%
%===================================================================================================

% Set default variable values as necessary
if isempty(cMapRange)
    cMapRange = [0 5000];
end
if isempty(sigma)
    sigma = 0.5;
end
if isempty(fps)
    fps = 30;
end

% Get movie data
disp('loading currData')
currData = sessionObj.wholeSession(:,:,planeNum,:,trialNum);
disp('currData loaded')
trialData = squeeze(currData);
filtData = [];
for iVol = 1:size(trialData,3)
    filtData(:,:,iVol) = imgaussfilt(trialData(:,:,iVol), sigma);
end

% Open movie in video player app
h = implay_scaled(filtData, cMapRange, fps);

% Adjust movie size and colormap scaling
h.Parent.Position = [700 300 820 600];
oldPos = h.Visual.Axes.Position;
drawnow
h.Visual.Axes.Position = ([0 1 3.2 3.2]  .*oldPos) - [0 0 0 0];
drawnow

end%function