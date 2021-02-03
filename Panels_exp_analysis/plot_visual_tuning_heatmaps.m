function f = plot_visual_tuning_heatmaps(sourceTbl, plotFl, plotParams)
%===================================================================================================
% Plots a heatmap for each EB wedge or PB glomerulus depicting mean visual tuning at each bar 
% location (x axis) during each trial (y axis). Also draws a box around trials in which a drug 
% was applied. SourceTbl needs to have exactly one row for each trial in plotFl.
%
% INPUTS:
%
%       sourceTbl   = source data table with one row for each trial to be plotted. Will use the 
%                     following specific fields: 
%                           trialNum
%                           startTime  (from drug timing metadata)
%                           visualStim (from drug timing metadata)
% 
%       plotFl      = an array of fluorescence data with dimensions [barPos, trial, EB wedge] 
%
%       plotParams  = struct containing the following fields of plotting parameters:
%                           smWin         (width of gaussian smoothing window for imaging data)
%                           flMax         (value to cap fl at for imagesc plots, or [] to skip)
%                           figPos        (width and height of the figure, or [] to skip)
%                           figNum        (number of figure to create, or [] to default to 1)
%                           SV, SH, ML, MR, MT, and MB (subaxis spacing arguments)
% 
% OUTPUTS:
%
%       f   = handle to the figure that was created
% 
%===================================================================================================

tbl = sourceTbl;
p = plotParams;

% Create figure
if isempty(p.figNum)
   p.figNum = 1;
end
f = figure(p.figNum);clf;
if ~isempty(p.figPos)
    f.Position(3:4) = p.figPos;
    if sum(f.Position([2, 4])) > 1080
        f.Position(2) = 50;
    end
end
f.Color = [1 1 1];

% Determine subplot grid size
nPlots = size(plotFl, 3);
if nPlots == 3
    plotPos = [2 2]; % Because [1 3] looks bad
else
    plotPos = numSubplots(nPlots);
end

drugTrials = find(~isnan(tbl.startTime));

for iPlot = 1:nPlots
    
    subaxis(plotPos(1), plotPos(2), iPlot, 'mt', p.MT, 'mb', p.MB, 'sv', p.SV, 'mr', p.MR, ...
            'ml', p.ML, 'sh', p.SH);
    
    currData = smoothdata(plotFl(:, :, iPlot), 1, 'gaussian', p.smWin);
    if ~isempty(p.flMax)
        currData(currData > p.flMax) = p.flMax;
    end
    
    % Plot data
    barAngles = (-180 + 3.75):3.75:180;
    xx = [barAngles(1), barAngles(end)];
    imagesc(xx, [1 size(currData, 2)], currData');
    hold on; plot([0 0], [0 size(currData, 2) + 1], 'color', 'k', 'linewidth', 2)
    ylim([0.5 size(currData, 2) + 0.5])
    xlabel('Bar position (degrees from front of fly)');
    ylabel('Trial');
    colormap('magma')
    
    % Add Y tick labels
    ax = gca;
    ax.YTick = 1:size(currData, 2);
    ax.YTickLabels = tbl.trialNum;
    
    % Remove extraneous X and Y axis labels
    if mod(iPlot - 1, plotPos(2))
        ax.YLabel.String = [];
    end
    if iPlot <= (nPlots - plotPos(2))
        ax.XLabel.String = [];
        ax.XTickLabel = [];
    end
    
    ax.FontSize = 12;
    
    % Indicate transitions to drug application trials
    if ~isempty(drugTrials)
        for iTrial = 1:numel(drugTrials)
            xL = xlim();
            if strcmp(tbl.visualStim{drugTrials(iTrial)}, 'yes')
                lineColor = 'g';
            else
                lineColor = rgb('yellow');
            end
            plot(xL, [-0.5, -0.5] + drugTrials(iTrial), 'color', lineColor, 'linewidth', 3)
            plot(xL, [0.5, 0.5] + drugTrials(iTrial), 'color', lineColor, 'linewidth', 3)
            plot([xL(1), xL(1)] + diff(xL)*0.005, [-0.5, +0.5] + drugTrials(iTrial), 'color', lineColor, ...
                    'linewidth', 3)
            plot([xL(2), xL(2)] - diff(xL)*0.005, [-0.5, +0.5] + drugTrials(iTrial), 'color', lineColor, ...
                    'linewidth', 3)
        end
    end
    
end%iPlot

end