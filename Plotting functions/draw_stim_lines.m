function hPlots = draw_stim_lines(stimEpochs, stimShadingColors, varargin)
%===================================================================================================

% For drawing vertical lines on a 2D-behavior annotation style plot to indicate stim times.

% INPUTS:
%       stimEpochs         = 
%
%       stimShadingColors  = Cell array containing strings specifying a line color for each row in 
%                            stimEpochs.
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       'plotAxes'    = (default: gca) 
%
%       'frameRate'   = (default: 25) 
%
%       'yLims'       = (default: ylim())
%
% OUTPUTS:
%       hPlots  = A cell array of the handles to each line that was drawn
% 
%===================================================================================================

% Parse optional arguments
p = inputParser;
addParameter(p, 'plotAxes', gca());
addParameter(p, 'frameRate', 25);
addParameter(p, 'yLims', []);
parse(p, varargin{:});
plotAxes = p.Results.plotAxes;
frameRate = p.Results.frameRate;
yLims = p.Results.yLims;
if isempty(yLims)
    yLims = ylim(plotAxes);
end

% Plot lines
hPlots = [];
nStimEpochs = size(stimEpochs, 1);
for iStim = 1:nStimEpochs
    stimOnsetFrame = stimEpochs(iStim, 1) * frameRate;
    stimOffsetFrame = stimEpochs(iStim, 2) * frameRate;
    shadeColor = rgb(stimShadingColors{iStim});
    hPlots{iStim, 1} = plot(plotAxes, [stimOnsetFrame, stimOnsetFrame], yLims, 'Color', shadeColor, 'LineWidth', 2);
    hPlots{iStim, 2} = plot(plotAxes, [stimOffsetFrame, stimOffsetFrame], yLims, 'Color', shadeColor, 'LineWidth', 2);
end





end