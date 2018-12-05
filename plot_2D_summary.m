function [plotHandle, plotAxes, plotFig] = plot_2D_summary(infoStruct, dataArr, varargin)
%===================================================================================================
% PLOT A 2D SUMMARY FIGURE OF DATA THROUGHOUT EXPERIMENT
%
% Uses imagesc() to visualize 2D data in a [trial, volume/frame] format from an entire experiment.
% By default will plot trials in chronological order, but they can be optionally split into groups.
% Can create new figure for the plot or plot in a specific axes object.
%
% INPUTS:
%       infoStruct    = the main imaging data structure containing metadata for the experiment.
%                         Specifically, must contain the fields "expDate", "trialDuration" and
%                         "volumeRate" (latter is only necessary if using default 'sampRate' value.
%
%       dataArr       = array of data to be plotted in format [trial, volume/frame]
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%       'plotAxes'    = (default: [])
%
%       'figPos'      = (default: [-1050 45 1050 950])
%
%       'trialGroups' = (default: [])
%
%       'titleStr'    = (default: infoStruct.expDate)
%
%       'saveDir'     = (default: 0)
%
%       'fileName'    = (default: 'Untitled')
%
%       'xAxLabel'    = (default: 'Time (sec)')
%
%       'colormap'    = (default: parula)
%
%       'sampRate'    = (default: infoStruct.volumeRate)
%
% OUTPUTS:
%       plotHandle  = the handle to the plot that was created by the function
%
%       plotAxes    = the axes that the figure was plotted in
%
%       plotfig     = handle to the new figure, if one was created (otherwise returns [])
%===================================================================================================

% Parse optional arguments
p = inputParser;
addParameter(p, 'plotAxes', []);
addParameter(p, 'figPos', [-1050 45 1050 950])
addParameter(p, 'trialGroups', [])
addParameter(p, 'titleStr', regexprep(infoStruct.expDate, '_', '\\_'));
addParameter(p, 'saveDir', 0);
addParameter(p, 'fileName', 'Untitled');
addParameter(p, 'xAxLabel', 'Time (sec)');
addParameter(p, 'colormap', parula(numel(unique(dataArr))));
addParameter(p, 'sampRate', infoStruct.volumeRate);
parse(p, varargin{:});
plotAxes = p.Results.plotAxes;
figPos = p.Results.figPos;
trialGroups = p.Results.trialGroups;
titleStr = p.Results.titleStr;
saveDir = p.Results.saveDir;
fileName = p.Results.fileName;
xAxLabel = p.Results.xAxLabel;
cm = p.Results.colormap;
sampRate = p.Results.sampRate;

% Create or select figure and axes depending on whether an axes handle was provided
if isempty(plotAxes)
    plotFig = figure(1); clf;
    plotAxes = axes();
    plotFig.Position = [-1050 45 1050 950];
    plotFig.Color = [1 1 1];
else
    axes(plotAxes)
    plotFig = gcf();
end

% Plot data
if isempty(trialGroups)
    plotHandle = imagesc(plotAxes, dataArr);
else
    % Separate into trial groups, adding black spacers in between
    plotArr = [];
    spacerArr = ones(4, size(dataArr, 2)) * (min(dataArr(:)) - 1);
    for iGroup = 1:numel(unique(trialGroups))
        if iGroup == 1
            plotArr = dataArr(trialGroups == iGroup, :);
        else
            plotArr = [plotArr; spacerArr; dataArr(trialGroups == iGroup, :)];
        end
    end
    plotHandle = imagesc(plotAxes, plotArr);
    
    % Update colormap to distinguish spacer rows
    colormap(plotAxes, [0 0 0; cm(1:end-1,:)])
end
colormap(plotAxes, cm)

% Format axes
ax = gca();
ax.FontSize = 14;
ax.YTick = 0:10:sum(infoStruct.goodTrials);
ax.XTick = [0:5:infoStruct.trialDuration] * sampRate;
ax.XTickLabel = 0:5:infoStruct.trialDuration;
title(titleStr);
xlabel(xAxLabel)
ylabel('Trial')

% Save figure
if saveDir
    export_fig(fullfile(saveDir, fileName), '-png', plotFig);
    if ~isdir(fullfile(saveDir, 'figFiles'))
        mkdir(fullfile(saveDir, 'figFiles'))
    end
    savefig(plotFig, fullfile(saveDir, 'figFiles', fileName));
    
end
end