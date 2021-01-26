function [f, ax] = plot_single_trial_bump_heatmaps(dataTbl, plotParams)
%===================================================================================================
% Plots a column of heatmaps of bump fluorescence over time for each trial in dataTbl.
%
% INPUTS:
%
%       dataTbl     = source data table with one row for each trial to be plotted. Will use the 
%                     following specific fields: 
%                               rawFl, trialDff, or expDff (depending on p.flType)
%                               expID
%                               trialNum
%                               trialDuration
%                               volTimes
%                               cycStartTimes
%                               drugName (this and subsequent fields are from drug timing metadata)
%                               startTime
%                               duration
%                               visualStim
% 
%       plotParams  = struct containing the following fields of plotting parameters:
%                           smWin         (width of gaussian smoothing window for imaging data)
%                           flType        ('rawFl', 'trialDff', or 'expDff')
%                           plotBarCycles (boolean, include vertical lines at bar cycle bounds?)
%                           flMax         (value to cap fl at for imagesc plots, or [] to skip)
%                           figPos        (width and height of the figure, or [] to skip)
%                           figNum        (number of figure to create, or [] to default to 1)
%                           SV, SH, ML, MR, MT, and MB (subaxis spacing arguments)
% 
% OUTPUTS:
%
%       f   = handle to the figure that was created
% 
%       ax  = array of axes handles for each subplot
% 
%===================================================================================================

p = plotParams;
tbl = dataTbl; 

% Get correct source data
sourceData = tbl.(p.flType);
smFl = cellfun(@(x) smoothdata(x, 1, 'gaussian', p.smWin), sourceData, 'uniformoutput', 0);

% Find global max and min fl values to keep color scale consistent
flRanges = cell2mat(cellfun(@(x) [min(x(:), [], 'omitnan'), max(x(:), [], 'omitnan')], smFl, ...
        'uniformoutput', 0));
plotFlRange = [min(flRanges(:), [], 'omitnan'), max(flRanges(:), [], 'omitnan')];

% Create and position figure
if isempty(p.figNum)
    p.figNum = 1;
end
f = figure(p.figNum); clf;
f.Color = [1 1 1];
if ~isempty(p.figPos)
    f.Position(3:4) = figPos;
    if sum(f.Position([2, 4])) > 1080
        f.Position(2) = 50;
    end
end
colormap(magma);

% Identify any drug treatment trials
drugTrials = find(~isnan(tbl.startTime));

for iTrial = 1:size(tbl, 1)
    
    % Plot fl data
    ax(iTrial) = subaxis(size(tbl, 1), 1, iTrial, 'mt', p.MT, 'mb', p.MB, 'sv', p.SV, 'mr', p.MR, ...
            'ml', p.ML, 'sh', p.SH);
    
    plotFl = smFl{iTrial};
    plotFl(1, 1) = plotFlRange(1);
    plotFl(end, end) = plotFlRange(2);
    if ~isempty(p.flMax)
        plotFl(plotFl > p.flMax) = p.flMax;
    end
    imagesc([0, tbl.trialDuration(iTrial)], [1, 8], plotFl');
    hold on
    
    % Plot PVA if using wedge data
    if size(plotFl, 2) == 8
        plot(tbl.volTimes{iTrial}, 9 - smoothdata(tbl.pvaWedge{iTrial}, 1, 'gaussian', p.smWin), ...
            'color', 'c', 'linewidth', 0.5)
    end
    
    % Plot any drug applications in the current trial
    if ~isempty(drugTrials) && ismember(iTrial, drugTrials)
        if strcmp(tbl.visualStim{iTrial}, 'yes')
            lineColor = 'g'; % Use green for trials with visual stim, yellow for darkness trials
        else
            lineColor = rgb('yellow');
        end
        yL = ylim();
        xx_1 = [1, 1] .* tbl.startTime(iTrial);
        xx_2 = [1, 1] .* tbl.startTime(iTrial) + tbl.duration(iTrial);
        plot(xx_1, yL, 'color', lineColor, 'linewidth', 5);
        plot(xx_2, yL, 'color', lineColor, 'linewidth', 5);
        ylim(yL)
    end
    
    % Plot bar cycle boundaries if necessary
    currStartVols = tbl.cycStartVols{iTrial};
    if p.plotBarCycles && ~isempty(currStartVols)
        yL = ylim();
        for iCycle = 1:numel(currStartVols)
            xx = [1, 1] .* tbl.volTimes{iTrial}(currStartVols(iCycle));
            plot(xx, yL, 'color', 'r', 'linewidth', 1)
        end
        ylim(yL);
    end
    
    % Label each plot with trial number
    ax(iTrial).YTickLabel = [];
    t = ylabel(['Trial #', num2str(tbl.trialNum(iTrial))], 'rotation', 0, 'FontSize', 12);
    t.HorizontalAlignment = 'right';
    t.VerticalAlignment = 'middle';
    t.Position(1) = t.Position(1) * 2;
    
    % Label X axis on final trial
    if iTrial < size(tbl, 1)
        ax(iTrial).XTickLabel = [];
    else
        str = ax(iTrial).YLabel.String;
        ax(iTrial).FontSize = 12;
        xlabel('Time (sec)', 'FontSize', 14);
        ylabel(str, 'rotation', 0, 'FontSize', 12);
    end
    
    % Add title above the first plot
    if iTrial == 1
        titleStr = [tbl.expID{iTrial}, '  —  single trial ', p.flType];
        if ~isempty(drugTrials)
            titleStr = [titleStr, '  (green/yellow lines = ', tbl.drugName{drugTrials(1)}, ...
                ' onset and offset)'];
        end
        title(titleStr, 'FontSize', 16);
    end   
    
end%iTrial





end%function