function f = plot_visual_tuning_curves_polar(sourceTbl, plotFl, plotParams)
%===================================================================================================
% Plots a tuning curves for each EB wedge or PB glomerulus overlaid on a polar plot for each trial. 
% Also indicates in the titles of the plots whether a drug was applied. SourceTbl needs to have 
% exactly one row for each trial in plotFl.
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
%                           matchRLims
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

drugTrials = find(~isnan(tbl.startTime));

nPlots = size(plotFl, 2);
for iPlot = 1:nPlots
    
    
    
    currData = smoothdata(squeeze(plotFl(:, iPlot, :)), 1, 'gaussian', p.smWin);
    nRois = size(currData, 2);
    
    % Create axes
    ax = subaxis(1, nPlots, iPlot, 'mt', p.MT, 'mb', p.MB, 'sv', p.SV, 'mr', p.MR, ...
        'ml', p.ML, 'sh', p.SH);
    hold on;
    pax = polaraxes;
    hold on;
    pax.Position = ax.Position;
    ax.Visible = 'off';
    
    % Plot data
    barAngles = deg2rad((-180 + 3.75):3.75:180);
        cm = hsv(nRois) .* 0.85;
%         cm = phasemap(nRois);
%     cm = colorcet('C2', 'N', nRois);
    for iWedge = 1:nRois
        %         plotData = currData(:, iWedge) - min(currData(:));
        plotData = currData(:, iWedge);
        polarplot(barAngles, plotData, 'color', cm(iWedge, :), 'linewidth', 1.5);
    end
    
   % Format axes
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.ThetaTick = [0:45:179, 225:45:359];
    pax.ThetaTickLabel = {'0' '+45', '+90' '+135' '-135', '-90', '-45'};
    pax.RTick = [];
    pax.FontSize = 12;
    if p.matchRLims
        pax.RLim(2) = 1.01 * max(as_vector(smoothdata(plotFl, 1, 'gaussian', p.smWin)));
    end
    title(['#', num2str(iPlot)], 'fontSize', 14)
    
    % Indicate transitions to drug application trials
    if ~isempty(drugTrials) && ismember(iPlot, drugTrials)
        pax.Title.String = [pax.Title.String, ' (post-', tbl.drugName{iPlot}, ')'];
    end
    
end%iPlot


end