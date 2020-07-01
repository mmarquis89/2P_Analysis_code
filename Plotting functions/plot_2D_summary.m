function [plotHandle, plotAxes, plotFig] = plot_2D_summary(dataArr, sampRate, varargin)
%===================================================================================================
% PLOT A 2D SUMMARY FIGURE OF DATA THROUGHOUT EXPERIMENT
%
% Uses imagesc() to visualize 2D data in a [trial, volume/frame] format from an entire experiment.
% By default will plot trials in chronological order, but they can be optionally split into groups.
% Can create new figure for the plot or plot in a specific axes object. Any trials for which 
% trialGroups == 0 will be omitted from the plots.
%
% INPUTS:
%
%       dataArr       = array of data to be plotted in format [trial, volume/frame]
%
%       sampRate      = the sampling rate of dataArr in samples/sec
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%       'plotAxes'      = (default: [])
%
%       'figPos'        = (default: [-1050 45 1050 950])
%
%       'trialGroups'   = (default: [])
%
%       'titleStr'      = (default: '')
%
%       'saveDir'       = (default: 0)
%
%       'fileName'      = (default: 'Untitled')
%
%       'xAxLabel'      = (default: 'Time (sec)')
%
%       'colormap'      = (default: parula)
%
%       'spacerArrRows' = (default: 4)
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
addParameter(p, 'titleStr', '');
addParameter(p, 'saveDir', 0);
addParameter(p, 'fileName', 'Untitled');
addParameter(p, 'xAxLabel', 'Time (sec)');
addParameter(p, 'colormap', parula(numel(unique(dataArr))));
addParameter(p, 'spacerArrRows', 4);
parse(p, varargin{:});
plotAxes = p.Results.plotAxes;
figPos = p.Results.figPos;
trialGroups = p.Results.trialGroups;
titleStr = p.Results.titleStr;
saveDir = p.Results.saveDir;
fileName = p.Results.fileName;
xAxLabel = p.Results.xAxLabel;
cm = p.Results.colormap;
spacerArrRows = p.Results.spacerArrRows;

% Create or select figure and axes depending on whether an axes handle was provided
if isempty(plotAxes)
    plotFig = figure(1); clf;
    plotAxes = axes();
    plotFig.Position = figPos;
    plotFig.Color = [1 1 1];
else
    axes(plotAxes)
    plotFig = gcf();
end

% Plot data
if isempty(trialGroups) || numel(unique(trialGroups(trialGroups > 0))) == 1
    if isempty(trialGroups) || numel(unique(trialGroups)) == 1
        plotHandle = imagesc(plotAxes, dataArr);
    else
        plotHandle = imagesc(plotAxes, dataArr(trialGroups > 0, :));
    end
    colormap(plotAxes, cm)
else
    % Separate into trial groups, adding black spacers in between
    plotArr = []; spacerInds = [];
    spacerArr = ones(spacerArrRows, size(dataArr, 2)) * max(dataArr(:)) * 1.01;
    for iGroup = 1:numel(unique(trialGroups(trialGroups > 0)))
        if iGroup == 1
            plotArr = dataArr(trialGroups == iGroup, :);
            spacerCount = zeros(1, size(dataArr(trialGroups == iGroup, :), 1));
        else
            spacerInds = [spacerInds, (size(plotArr, 1)+1):(size(plotArr, 1)+size(spacerArr, 1))];
            plotArr = [plotArr; spacerArr; dataArr(trialGroups == iGroup, :)];
            spacerCount = [spacerCount, ones(1, size(spacerArr, 1) + ...
                    size(dataArr(trialGroups == iGroup, :), 1)) * (iGroup - 1)];
        end
        
    end
    plotHandle = imagesc(plotAxes, plotArr);
    
    % Update colormap to distinguish spacer rows
    colormap(plotAxes, [cm(1:end,:); 0 0 0])
end

yTickNums = 0:10:size(dataArr, 1);
yTickLocs = yTickNums;
if ~isempty(trialGroups) && numel(unique(trialGroups(trialGroups > 0))) > 1
    for iTick = 1:numel(yTickNums)
        if yTickLocs(iTick) > 0
           offset =  spacerCount(yTickLocs(iTick)) * spacerArrRows;
        else
            offset = 0;
        end
        yTickLocs(iTick) = yTickLocs(iTick) + offset;
    end
end

% Format axes
ax = gca();
ax.FontSize = 12;
ax.YTick = yTickLocs;
ax.YTickLabel = yTickNums;
trialDuration = size(dataArr, 2) / sampRate;
ax.XTick = (0:5:trialDuration) * sampRate;
ax.XTickLabel = 0:5:trialDuration;
title(titleStr);
xlabel(xAxLabel)
ylabel('Trial')

% % Update Y ticks to account for spacers and omitted trials
% if ~isempty(trialGroups) %&& numel(unique(trialGroups(trialGroups > 0))) > 1
%     testTicks = [];
%     for iTick = 1:numel(ax.YTick)
%         testTicks(iTick) = ax.YTick(iTick) + ...
%             (sum(spacerInds(1:4:end) < ax.YTick(iTick))*4) - ...
%             sum(~trialGroups(1:ax.YTick(iTick)));
%     end
%     ax.YTick = testTicks;
%     ax.YTickLabel = 0:10:infoStruct.nTrials;%(testTicks(end) + sum(~trialGroups));
% end

% Save figure
if saveDir
    export_fig(fullfile(saveDir, fileName), '-png', plotFig);
    if ~isfolder(fullfile(saveDir, 'figFiles'))
        mkdir(fullfile(saveDir, 'figFiles'))
    end
    savefig(plotFig, fullfile(saveDir, 'figFiles', fileName));
    
end
end