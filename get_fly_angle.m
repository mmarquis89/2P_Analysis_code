function [thetaDeg, thetaRad] = get_fly_angle()
%=======================================================================================================
% MEASURES ANGLE OF FLY ON FICTRAC BALL
%
% Allows user to select an image of a fly on a FicTrac ball viewed directly from the side and measure
% the angle that the fly is offset from standing directly on top of the ball. This angle can then be
% used to manually edit the fictrac transform matrix to ensure accurate tracking.
%
% OUTPUTS:
%       
%       thetaDeg = fly's offset angle in degrees (for easy sanity checking)
%
%       thetaRad = the same angle in radians (this value should be used in the transform file)
%
%========================================================================================================

% Get and load image file
[fileName, pathName] = uigetfile('*');
if fileName == 0
    disp('Angle selection cancelled')
    thetaDeg = [];
    thetaRad = [];
    return
end

% Read image (or extract first frame if image is a video)
if regexp(fileName, '.*(.avi|.mp4|.mov)')
    myVid = VideoReader(fullfile(pathName, fileName));
    firstFrame = readFrame(myVid);
    img = uint8(firstFrame(:,:,1));
else
    img = uint8(imread(fullfile(pathName, fileName)));
end

% Create figure and plot image
f = figure(1); clf
imshow(img, []);
ax = gca();
hold on;

% Get 3 points on circumference of circle from user
ax.Title.String = 'Select 3 points on the circumference of the ball'
[xC, yC] = ginput(3);

% Plot circle
ABC = [xC, yC];
[R, xcyc] = fit_circle_through_3_points(ABC);
xStart = xcyc(1);
yStart = xcyc(2);

% Get point directly under fly from user
ax.Title.String = 'Select a point on the ball directly below the fly'
[flyX, flyY] = ginput(1);

% Calculate angle from center of circle
flyVec = [flyX, flyY] - [xStart, yStart];
[thetaRad, ~] = cart2pol(flyVec(1), flyVec(2));
thetaRad = thetaRad + pi/2; % rotating 90° so zero is up 
thetaDeg = rad2deg(thetaRad);

close(f);

end%function






