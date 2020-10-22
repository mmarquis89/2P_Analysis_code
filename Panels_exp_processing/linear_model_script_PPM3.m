%% Set up a new model object

% Set source data parameters;
p = [];
p.maxSpeed = 100;
p.smWinVols = 3;
p.roiName = 'EB-DAN';
p.smWinFrames = 5;
p.smReps = 20;
p.ftLagVols = 3;
p.speedType = 'fwSpeed';

expIDList = {'20201015-1', '20201015-2', '20201019-1', '20201019-2'};

skipTrials = {[11:15], [11:16], [8:15], [1:4, 9, 12:13]};

skipVols = repmat({[]}, 1, numel(skipTrials));
                     
         
expInfoTbl = table(expIDList', skipTrials', skipVols', 'VariableNames', {'expID', 'skipTrials', ...
        'skipVols'});

% Create analysis object
rm = RegressionModelAnalysis_PPM3(expInfoTbl, p);

%% Choose parameters and initialize model
mp = [];
mp.trainTestSplit = 0.8;
mp.criterion = 'adjrsquared'; % 'sse, 'aic', 'bic', '', or 'adjrsquared'
mp.upper = ['linear'];
mp.pEnter = [0.03];
mp.pRemove = [0];
mp.verbose = 0;
mp.useYaw = 1;
mp.useDriftCorrection = 0;

% Initialize models
rm = rm.initialize_models(mp);

%% Train initial models
try
fullMdls = {};
fullMdlPredFl = {};
fullMdlAdjR2 = [];
disp('Training initial models...')
for iExp = 1:size(rm.sourceData, 1)
    if numel(rm.modelParams) > 1
        mp = rm.modelParams(iExp);
    else
        mp = rm.modelParams;
    end
    disp(rm.sourceData.expID{iExp})
    
    % Select and split up data from current experiment
    currModelData = rm.modelData(iExp, :);
    tblFit = currModelData.fullDataTbl{:}(currModelData.fitRowInds{:}, :);
    tblTest = currModelData.fullDataTbl{:}(currModelData.testRowInds{:}, :);
    tblPred = currModelData.fullDataTbl{:}; 
    
    % Use all training data to create and evaluate a stepwise model
    kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
            mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
    emptyArgs = cellfun(@isempty, kvArgs);
    kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
    fullMdls{iExp} = stepwiselm(tblFit, kvArgs{:});

    [predFl, ~] = predict(fullMdls{iExp}, tblTest(~logical(sum(isnan(table2array(tblTest)), 2)), ...
            1:end-1));
    [~, fullMdlAdjR2(iExp)] = r_squared(tblTest.fl(~logical(sum(isnan(table2array(tblTest)), 2))), ...
            predFl, fullMdls{iExp}.NumCoefficients);
    [fullMdlPredFl{iExp}, ~] = predict(fullMdls{iExp}, ...
            tblPred(~logical(sum(isnan(table2array(tblPred)), 2)), 1:end-1));
end
disp('Final models created');
rm.modelData.fullMdls = fullMdls';
rm.modelData.fullMdlPredFl = fullMdlPredFl';
rm.modelData.fullMdlAdjR2 = fullMdlAdjR2';
disp(rm.modelData)
catch ME; rethrow(ME); end

%% Plot summary grid of coefficients used in the models

try
allCoeffNames = {};
for iExp = 1:size(rm.modelData, 1)
    allCoeffNames = [allCoeffNames, rm.modelData.fullMdls{iExp}.CoefficientNames];    
end
uniqueCoeffs = unique(allCoeffNames);
uniqueCoeffs = uniqueCoeffs(2:end); % Drop intercept term
%%
saveFig = 0;
uniqueCoeffs = uniqueCoeffs;%([1 2]);
coeffArrDc = zeros(numel(uniqueCoeffs), size(rm.modelData, 1));
for iExp = 1:size(rm.modelData, 1) 
    for iCoeff = 1:numel(uniqueCoeffs)
        if ismember(uniqueCoeffs{iCoeff}, rm.modelData.fullMdls{iExp}.CoefficientNames)
            coeffArrDc(iCoeff, iExp) = rm.modelData.fullMdls{iExp}.Coefficients.Estimate(...
                    strcmp(rm.modelData.fullMdls{iExp}.CoefficientNames, ...
                    uniqueCoeffs{iCoeff}));
        end
    end
end
coeffTblDc = array2table(coeffArrDc', 'variableNames', uniqueCoeffs, 'rowNames', ...
        rm.modelData.expID);
f = figure(37); clf; 
f.Color = [1 1 1];
f.Position = [f.Position(1:2), 835,950];
if (f.Position(2) + f.Position(4)) > 1000
    f.Position(2) = 1000 - f.Position(4);
end
plotArr = table2array(coeffTblDc);
lim = max(abs([max(plotArr(:)), min(plotArr(:))]));
xSize = size(plotArr, 2);
plotArr(:, end + 1) = -lim;
plotArr(:, end + 1) = lim;
ax = subaxis(1, 1, 1, 'ml', 0.14, 'mt', 0.12, 'mb', 0.03, 'mr', 0.11);
imagesc(plotArr);
cm = colormap(bluewhitered);
xlim([0.5, xSize + 0.5])
ax.XTick = 1:numel(uniqueCoeffs);
ax.XTickLabel = uniqueCoeffs;
ax.XTickLabelRotation = 30;
ax.YTick = 1:size(coeffTblDc, 1);
% ax.YTickLabel = cellfun(@(x, y) [x, ' (R^2=', num2str(y, 2), ')'], rm.modelData.expID, ...
%         mat2cell(rm.modelData.driftCorrectedMdlAdjR2, ones(1, size(rm.modelData, 1))), ...
%         'uniformoutput', 0);
ax.YTickLabel = rm.modelData.expID;
ax.XAxisLocation = 'top';
ax.FontSize = 14;
pos = ax.Position;
cb = colorbar();
cb.Position(1) = sum(ax.Position([1 3])) + 0.03;
ax.Position = pos;
hold on; 
xL = xlim();
yL = ylim();
for iCol = 1:(numel(uniqueCoeffs))
    if iCol == 1
        plot([iCol iCol] - 0.5, yL, 'color', 'k', 'linewidth', 1.5)
    end
    plot([iCol iCol] + 0.5, yL, 'color', 'k', 'linewidth', 1.5)
end
for iRow = 1:(numel(rm.modelData.expID))
    if iRow == 1
        plot(xL, [iRow iRow] - 0.5, 'color', 'k', 'linewidth', 1.5)
    end
    plot(xL, [iRow iRow] + 0.5, 'color', 'k', 'linewidth', 1.5)
    text(1.5, iRow, ['   Adj. R^2 = ', num2str(rm.modelData.fullMdlAdjR2(iRow), 2)], ...
            'FontSize', 12);
end
title('Drift-corrected model coefficients')
f.UserData.sourceDataParams = rm.sourceDataParams;
f.UserData.modelParams = rm.modelParams;
if saveFig
    pEnterStr = num2str(rm.modelParams.pEnter);
    saveFileName = ['model_coeff_summary_moveSpeed_noYaw_withInteractions'];
    save_figure(f, figDir, saveFileName);
end
catch ME; rethrow(ME); end

%% Plot predicted vs. measured fluorescence

saveFigs = 0;
fileNameStr = 'fullExp_';

predictorVars = {'fwSpeed', 'yawSpeed'}; % odorResp, fwSpeed, moveSpeed, yawSpeed
legendLocs = {'sw', 'sw', 'nw', 'nw'};
% legendLocs = {'best', 'best', 'best', 'best'};
xLimits = repmat({[]}, 1, 4);

try
legendStr = {'predicted fl', 'measured fl'};
for iExp = 1:size(rm.modelData, 1)
    
    f = figure(iExp + 100); clf;
    f.Color = [1 1 1];
    
    currModelData = rm.modelData(iExp, :);
    currSourceData = rm.sourceData(iExp, :);
       
    fullDataTbl = currModelData.fullDataTbl{:};
    fullDataTbl = fullDataTbl(~isnan(fullDataTbl.fl), :);
    
    legendStr = {'predicted fl', 'measured fl'};
    predFl = currModelData.fullMdlPredFl{:};
    measuredFl = currSourceData.fl{:};
    xx = currSourceData.volTimes{:}(1:numel(measuredFl))';
    xL = xLimits{iExp};
    if isempty(xL)
        xL = [0 xx(end)];
    end
    
    % Measured Fl data
    ax = subaxis(24, 1, 7:16, 'mt', 0.05, 'mb', 0.07, 'ml', 0.03, 'mr', 0.03, 'sv', 0.01); hold on;
    plot(xx, measuredFl, 'linewidth', 1, 'color', rgb('green'));
    ax.YTickLabel = [];
    legend(legendStr{2}, 'fontsize', 12, 'Location', legendLocs{1}, 'autoupdate', 'off', 'box', ...
            'off');
    ax.FontSize = 14;
    ax.YTickLabel = [];
    ax.XTickLabel = [];
    ax.XLim = xL;
    ax.YLim = [min(measuredFl) max(measuredFl)];
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    box off;
    allAxes = ax;
    
    
    % Predicted Fl data
    ax = subaxis(24, 1, 17:24, 'mt', 0.05, 'mb', 0.07, 'ml', 0.03, 'mr', 0.03, 'sv', 0.01); hold on;
    plot(xx, predFl, 'linewidth', 1, 'color', rgb('orangered'));
    ax.YTickLabel = [];
    legend(legendStr{1}, 'fontsize', 12, 'Location', legendLocs{2}, 'autoupdate', 'off', 'box', ...
            'off');
    ax.FontSize = 14;
    ax.YTickLabel = [];
    ax.XTickLabel = [];
    ax.XLim = xL;
    ax.YLim = [min(predFl) max(predFl)];
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    box off;
    allAxes(end + 1) = ax;
    
    % Predictor variables
    ax = subaxis(24, 1, 1:6); hold on
    legendStr = {};
    if any(strcmpi('fwSpeed', predictorVars))
        plot(xx, fullDataTbl.fwSpeed(~logical(sum(isnan(table2array(fullDataTbl)), 2))), '-', ...
                'color', rgb('blue'), 'linewidth', 1);
        legendStr{end + 1} = 'abs(fwSpeed)';
    end
    if any(strcmpi('moveSpeed', predictorVars))
        plot(xx, fullDataTbl.moveSpeed(~logical(sum(isnan(table2array(fullDataTbl)), 2))), '-', ...
                'color', rgb('blue'), 'linewidth', 1);
        legendStr{end + 1} = 'moveSpeed)';
    end
    if any(strcmpi('yawSpeed', predictorVars))
        plot(xx, fullDataTbl.yawSpeed(~logical(sum(isnan(table2array(fullDataTbl)), 2))), 'color', ...
                rgb('red'));
        legendStr{end + 1} = 'yawSpeed';
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
        suptitle([currModelData.expID{:}, '  —  drift-corrected model (adj. R^2 = ', ...
            num2str(currModelData.fullMdlAdjR2, 2), ')'])
    allAxes(end + 1) = ax;
    
    linkaxes(allAxes, 'x')
    
    if saveFigs
        figTitle = ['predictedVsMeasuredFl_', fileNameStr, currModelData.expID{:}];
        save_figure(f, figDir, figTitle);
    end
end

catch ME; rethrow(ME); end