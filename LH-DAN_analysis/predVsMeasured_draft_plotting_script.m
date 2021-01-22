
%% Plot predicted vs. measured fluorescence
try
predictorVars = {'odorHistory'}; % odorHistory, odorResp, fwSpeed, yawSpeed
legendLocs = {'sw', 'sw', 'sw'};
xL = [];

legendStr = {'predicted fl', 'measured fl'};
for iExp = [1:4 7:size(rm.modelData, 1)]%;1:size(rm.modelData, 1)
    
    f = figure(iExp + 100); clf;
    f.Color = [1 1 1];
    
    currModelData = rm.modelData(iExp, :);
    currSourceData = rm.sourceData(iExp, :);
    
    bestWinSize = currModelData.bestWinSizes;
    
    fullDataTbl = currModelData.fullDataTbl{:};
    fullDataTbl = fullDataTbl(~isnan(fullDataTbl.fl), :);
    
    
    ax = subaxis(12, 1, 1:5, 'mt', 0.05, 'mb', 0.07, 'ml', 0.03, 'mr', 0.03, 'sv', 0.01); hold on;
    legendStr = {'measured fl', 'predicted fl (drift-corrected model)', };
    dcPredFl = smoothdata(currModelData.driftCorrectedMdlPredFl{:}, 'gaussian', 3);
    nohPredFl = currModelData.noOdorHistMdlPredFl{:};
    measuredFl = currModelData.driftCorrectedMdlMeasuredFl{:};
    xx = currSourceData.volTimes{:}(1:numel(measuredFl))';
    if isempty(xL)
        xL = [0 xx(end)];
    end
    plot(xx, measuredFl, 'linewidth', 1, 'color', rgb('green'));
    plot(xx, dcPredFl, 'linewidth', 1, 'color', rgb('orangered'));
    ax.YTickLabel = [];
    legend(legendStr, 'fontsize', 14, 'Location', legendLocs{1}, 'autoupdate', 'off', 'box', 'off');
    ax.FontSize = 14;
    ax.YTickLabel = [];
    ax.XTickLabel = [];
    ax.XLim = xL;
%     ax.YLim = [min(measuredFl) max(measuredFl)];
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    title([currModelData.expID{:}, '  —  drift-corrected model: adj. R^2 = ', ...
            num2str(currModelData.driftCorrectedMdlAdjR2, 2), '  (', num2str( ...
            currModelData.noOdorHistMdlAdjR2, 2), ' without odor history terms)'])
    box off;
    
    % Integrated odor history
    ax = subaxis(12, 1, 6:7); hold on
    legendStr = {};
    if any(strcmpi('odorResp', predictorVars))
        plot(xx, fullDataTbl.odorResp, 'color', rgb('red'), 'linewidth', 1);
        legendStr{end + 1} = 'Odor response';
    end
    if any(strcmpi('fwSpeed', predictorVars))
        plot(xx, fullDataTbl.fwSpeed, '-', 'color', rgb('blue'), 'linewidth', 1);
        legendStr{end + 1} = 'abs(fwSpeed)';
    end
    if any(strcmpi('yawSpeed', predictorVars))
        plot(xx, fullDataTbl.yawSpeed, 'color', rgb('red'));
        legendStr{end + 1} = 'yawSpeed';
    end
    if any(strcmpi('odorHistory', predictorVars))
        varName = ['odorHistory_', num2str(bestWinSize)];
        plot(xx, fullDataTbl.(varName), '-', 'color', rgb('black'), 'linewidth', 3);
        legendStr{end + 1} = regexprep(varName, '_', '\\_');
    end
    ax.FontSize = 14;
    ax.YTickLabel = [];
    ax.XTickLabel = [];
    ax.XLim = xL;
    box off
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    legend(legendStr, 'fontsize', 14, 'Location', legendLocs{2}, 'autoupdate', 'off', 'box', 'off');
    xlabel('Time (sec)')
    
    % No-odor history model
    ax = subaxis(12, 1, 8:12); hold on
    plot(xx, measuredFl, 'linewidth', 1, 'color', rgb('green'));
    plot(xx, nohPredFl, 'linewidth', 1, 'color', rgb('darkmagenta'));
    ax.FontSize = 14;
%     ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.XLim = xL;
    legendStr = {'measured fl', 'predicted fl (no odor history model)'};
    legend(legendStr, 'fontsize', 14, 'Location', legendLocs{3}, 'autoupdate', 'off', 'box', 'off');
    box off
%     ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
%     ax.YLim(2) = max(measuredFl);
    xlabel('Time (sec)')
    
end

catch ME; rethrow(ME); end

%% Plot predicted vs. measured fluorescence
try
predictorVars = {'odorHistory', 'odorResp', 'fwSpeed'}; % odorHistory, odorResp, fwSpeed, yawSpeed
legendLocs = {'sw', 'sw', 'w', 'w'};
xL = [];

legendStr = {'predicted fl', 'measured fl'};
for iExp = 6%1:size(rm.modelData, 1)
    
    f = figure(iExp + 10); clf;
    f.Color = [1 1 1];
    
    currModelData = rm.modelData(iExp, :);
    currSourceData = rm.sourceData(iExp, :);
    
    bestWinSize = currModelData.bestWinSizes;
    
    fullDataTbl = currModelData.fullDataTbl{:};
    fullDataTbl = fullDataTbl(~isnan(fullDataTbl.fl), :);
    
    legendStr = {'predicted fl (drift-corrected model)', 'predicted fl (no odor hist model)', ...
            'measured fl'};
    dcPredFl = currModelData.driftCorrectedMdlPredFl{:};
    nohPredFl = currModelData.noOdorHistMdlPredFl{:};
    measuredFl = currModelData.driftCorrectedMdlMeasuredFl{:};
    xx = currSourceData.volTimes{:}(1:numel(measuredFl))';
    if isempty(xL)
        xL = [0 xx(end)];
    end

    
    % Measured Fl data
    ax = subaxis(24, 1, 4:10, 'mt', 0.05, 'mb', 0.07, 'ml', 0.03, 'mr', 0.03, 'sv', 0.01); hold on;
    plot(xx, measuredFl, 'linewidth', 1, 'color', rgb('green'));
%     plot(xx, dcPredFl, 'linewidth', 1, 'color', rgb('orangered'));
    ax.YTickLabel = [];
    legend(legendStr{3}, 'fontsize', 12, 'Location', legendLocs{1}, 'autoupdate', 'off', 'box', ...
            'off');
    ax.FontSize = 14;
    ax.YTickLabel = [];
    ax.XTickLabel = [];
    ax.XLim = xL;
    ax.YLim = [min(measuredFl) max(measuredFl)];
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
%     title([currModelData.expID{:}, '  —  drift-corrected model (adj. R^2 = ', ...
%             num2str(currModelData.driftCorrectedMdlAdjR2, 2), ')'])
    box off;
    
    % Drift-corrected model
    ax = subaxis(24, 1, 11:17, 'mt', 0.05, 'mb', 0.07, 'ml', 0.03, 'mr', 0.03, 'sv', 0.01); hold on;
%     plot(xx, measuredFl, 'linewidth', 1, 'color', rgb('green'));
    plot(xx, dcPredFl, 'linewidth', 1, 'color', rgb('orangered'));
    ax.YTickLabel = [];
    legend(legendStr{1}, 'fontsize', 12, 'Location', legendLocs{2}, 'autoupdate', 'off', 'box', ...
            'off');
    ax.FontSize = 14;
    ax.YTickLabel = [];
    ax.XTickLabel = [];
    ax.XLim = xL;
    ax.YLim = [min(measuredFl) max(measuredFl)];
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
%     title([currModelData.expID{:}, '  —  drift-corrected model (adj. R^2 = ', ...
%             num2str(currModelData.driftCorrectedMdlAdjR2, 2), ')'])
    box off;
    
    % No-odor history model
    ax = subaxis(24, 1, 18:24); hold on
%     plot(xx, measuredFl, 'linewidth', 1, 'color', rgb('green'));
    plot(xx, nohPredFl, 'linewidth', 1, 'color', rgb('darkmagenta'));
    ax.FontSize = 14;
%     ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.XLim = xL;
    legend(legendStr{2}, 'fontsize', 12, 'Location', legendLocs{3}, 'autoupdate', 'off', ...
            'box', 'off');
    box off
%     ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    ax.YLim = [min(measuredFl) max(measuredFl)];
    
    % Integrated odor history
    ax = subaxis(24, 1, 1:3); hold on
    legendStr = {};
    if any(strcmpi('odorResp', predictorVars))
        plot(xx, fullDataTbl.odorResp, 'color', rgb('red'), 'linewidth', 1);
        legendStr{end + 1} = 'Odor response';
    end
    if any(strcmpi('fwSpeed', predictorVars))
        plot(xx, fullDataTbl.fwSpeed, '-', 'color', rgb('blue'), 'linewidth', 1);
        legendStr{end + 1} = 'abs(fwSpeed)';
    end
    if any(strcmpi('yawSpeed', predictorVars))
        plot(xx, fullDataTbl.yawSpeed, 'color', rgb('red'));
        legendStr{end + 1} = 'yawSpeed';
    end
    if any(strcmpi('odorHistory', predictorVars))
        varName = ['odorHistory_', num2str(bestWinSize)];
        plot(xx, fullDataTbl.(varName), '-', 'color', rgb('black'), 'linewidth', 3);
        legendStr{end + 1} = regexprep(varName, '_', '\\_');
    end
    ax.FontSize = 14;
    ax.YTickLabel = [];
    ax.XTickLabel = [];
    ax.XLim = xL;
    box off
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    legend(legendStr, 'fontsize', 12, 'Location', legendLocs{4}, 'autoupdate', 'off', 'box', ...
            'off');
    xlabel('Time (sec)')
    
end

catch ME; rethrow(ME); end