function h = implay_scaled(vidArr, mapRange, fps)
%===================================================================================================
% Calls the built-in implay() function, but adds the ability for the user to specify the range of 
% the colormap so it does not have to be changed manually.
%
% INPUTS:
%
%   vidArr    = an array containing the video data (M x 3), [y, x, frame]. You can also pass a 4-D
%               array as long as one dimension is a singleton, which will be squeezed out.
%
%   mapRange  = a 1x2 vector with the [min, max] values for the colormap scaling. Can pass a scalar 
%               value instead to use as the max map value and set the min value to zero.
%
%   fps       = the frame rate to play the video at.
%
% OUTPUTS:
%
%   h         = the handle to the video player app
%===================================================================================================

% Parse colormap argument
if length(mapRange) == 1
   mapRange = [0 mapRange];
elseif length(mapRange) ~= 2
    error('Error: invalid colormap range values')
end

% Open video
h = implay(squeeze(vidArr), fps);

% Apply range parameters
h.Visual.ColorMap.UserRangeMin = mapRange(1);
h.Visual.ColorMap.UserRangeMax = mapRange(2);
h.Visual.ColorMap.UserRange = 1;

end
%%%

