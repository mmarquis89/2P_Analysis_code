function plot_stim_shading(shadeTimes, varargin)
%===================================================================================================
% SHADE PLOT DURING STIMULUS PRESENTATION
%
% Adds a colored shading between a certain X interval on a plot, useful for indicating stimulus 
% presentations but could have other uses as well.
% 
% INPUTS:
%
%       shadeTimes = two element vector specifying the [onset offset] X-values of the shading
%
%
% OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%
%       'Axes'      = (default: gca) the handle of the axes to plot the shading in
% 
%       'Color'     = (default: [1 0 0]) rgb value for the shading color
% 
%       'Alpha'     = (default: 0.1) transparency value for the shading (ranging from 0-1).
%
%===================================================================================================

% Parse optional arguments
p = inputParser;
addParameter(p, 'Axes', gca());
addParameter(p, 'Color', [1 0 0]);
addParameter(p, 'Alpha', 0.1)
parse(p, varargin{:});
ax = p.Results.Axes;
fillColor = p.Results.Color;
fillAlpha = p.Results.Alpha;
onsetTime = shadeTimes(1);
offsetTime = shadeTimes(2);

% Huge Y values are so the bars don't get cut off if I increase the ylims later
yL = ylim(ax);
fill(ax, [onsetTime, onsetTime, offsetTime, offsetTime], [-1000, 1000, 1000, -1000], fillColor, ...
     'facealpha', fillAlpha, 'edgealpha', 0); 
ylim(ax, yL);

end