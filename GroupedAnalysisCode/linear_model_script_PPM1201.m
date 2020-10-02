parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
dataDir = fullfile(parentDir, 'all_experiments');
figDir = ...
        'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Figs\PPM1201_regression_analysis';

%% Set up a new model object

% Set source data parameters;
p = [];
p.maxSpeed = 100;
p.smWinVols = 5;
p.roiName = 'TypeF';
p.smWinFrames = 3;
p.smReps = 10;
p.ftLagVols = 3;
p.speedType = 'moveSpeed';
p.odorRespFilterDur = [7 5];

% Type F
expIDList = {'20180329-2', '20180405-2', '20180414-1', '20180414-2', '20180416-1', '20180523-2', ...
        '20181020-1', '20190226-3'};

skipTrials = {[], [], [], [], ...
              [], [], [], []};

skipVols = repmat({[]}, 1, numel(skipTrials));
                     
         
expInfoTbl = table(expIDList', skipTrials', skipVols', 'VariableNames', {'expID', 'skipTrials', ...
        'skipVols'});

% Create analysis object
rm = RegressionModelAnalysis_PPM1201(expInfoTbl, p);

%% Choose parameters and initialize model
mp = [];
mp.trainTestSplit = 0.8;
mp.criterion = 'adjrsquared'; % 'sse, 'aic', 'bic', '', or 'adjrsquared'
mp.upper = [];
mp.pEnter = [0.03];
mp.pRemove = [0];
mp.verbose = 0;
mp.useYaw = 0;
mp.useDriftCorrection = 1;

% Initialize models
rm = rm.initialize_models(mp);


%% =================================================================================================
%  MODEL PROCESSING
%  =================================================================================================

%% Train initial models
try
fullMdls = {};
fullMdlPredFl = {};
fullMdlAdjR2 = [];
disp('Training innitial models...')
for iExp = 1:size(rm.sourceData, 1)
    if numel(rm.modelParams) > 1
        mp = rm.modelParams(iExp);
    else
        mp = rm.modelParams;
    end
    disp(rm.sourceData.expID{iExp})
    
    % Select and split up data from current experiment
    currModelData = rm.modelData(iExp, :);
%     currModelData.fullDataTbl{:} = currModelData.fullDataTbl{:}(:, [1 2 4]);
    tblFit = currModelData.fullDataTbl{:}(currModelData.fitRowInds{:}, :);
    tblTest = currModelData.fullDataTbl{:}(currModelData.testRowInds{:}, :);
    tblPred = currModelData.fullDataTbl{:}; 
    
    % Use all training data to create and evaluate a stepwise model
    kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
            mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
    emptyArgs = cellfun(@isempty, kvArgs);
    kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
    fullMdls{iExp} = stepwiselm(tblFit, kvArgs{:});
    
    fullMdls{iExp} = fitlm(tblFit, 'fl ~ 1 + fwSpeed + odorResp + volsFromExpStart + fwSpeed:odorResp');

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

% Subtract volsFromOdorStart from fl data and re-fit models to "drift-corrected" data
try
disp('Training drift-corrected models...')
driftCorrectedMdls = {};
driftCorrectedMdlMeasuredFl = {};
driftCorrectedMdlTestFl = {};
driftCorrectedMdlPredFl = {};
driftCorrectedMdlAdjR2 = [];

for iExp = 1:size(rm.sourceData, 1)
    if numel(rm.modelParams) > 1
        mp = rm.modelParams(iExp);
    else
        mp = rm.modelParams;
    end
    disp(rm.sourceData.expID{iExp})
    
    % Select data from the current experiment
    currModelData = rm.modelData(iExp, :);
    tblPred = currModelData.fullDataTbl{:};
        
    currMdl = rm.modelData.fullMdls{iExp};
    if ismember('volsFromExpStart', currMdl.Coefficients.Properties.RowNames)
        tblPredVolAdjust = tblPred;
        tblPredVolAdjust.fl = tblPredVolAdjust.fl - (currMdl.Coefficients{'volsFromExpStart', ...
            'Estimate'} .* tblPredVolAdjust.volsFromExpStart);
        tblPred = tblPredVolAdjust(:, ~strcmp(tblPredVolAdjust.Properties.VariableNames, ...
                'volsFromExpStart'));
    end      
    tblFit = tblPred(currModelData.fitRowInds{:}, :);
    tblTest = tblPred(currModelData.testRowInds{:}, :);
    driftCorrectedMdlMeasuredFl{iExp} = tblPred.fl(~logical(sum(isnan(table2array(tblPred)), 2)));
    driftCorrectedMdlTestFl{iExp} = tblTest.fl(~logical(sum(isnan(table2array(tblTest)), 2)));
    
    % Use all training data to create and evaluate a final stepwise model
    kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
        mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
    emptyArgs = cellfun(@isempty, kvArgs);
    kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
    driftCorrectedMdls{iExp} = stepwiselm(tblFit, kvArgs{:});
    driftCorrectedMdls{iExp} = fitlm(tblFit, 'fl ~ 1 + fwSpeed + odorResp + fwSpeed:odorResp');

    [predFl, ~] = predict(driftCorrectedMdls{iExp}, ...
            tblTest(~logical(sum(isnan(table2array(tblTest)), 2)), 1:end-1));
    [~, driftCorrectedMdlAdjR2(iExp)] = r_squared(driftCorrectedMdlTestFl{iExp}, ...
            predFl, driftCorrectedMdls{iExp}.NumCoefficients);
    [driftCorrectedMdlPredFl{iExp}, ~] = predict(driftCorrectedMdls{iExp}, ...
            tblPred(~logical(sum(isnan(table2array(tblPred)), 2)), 1:end-1));
end
disp('Final models created');
rm.modelData.driftCorrectedMdls = driftCorrectedMdls';
rm.modelData.driftCorrectedMdlMeasuredFl = driftCorrectedMdlMeasuredFl'; 
rm.modelData.driftCorrectedMdlTestFl = driftCorrectedMdlTestFl';
rm.modelData.driftCorrectedMdlPredFl = driftCorrectedMdlPredFl';
rm.modelData.driftCorrectedMdlAdjR2 = driftCorrectedMdlAdjR2';
disp(rm.modelData)

catch ME; rethrow(ME); end

%% TEMPORARY HACK to remove volsFromExpStart term from all fullModels
fullMdls = {};
fullMdlPredFl = {};
fullMdlAdjR2 = [];
disp('Training innitial models...')
for iExp = 1:size(rm.sourceData, 1)
    if numel(rm.modelParams) > 1
        mp = rm.modelParams(iExp);
    else
        mp = rm.modelParams;
    end
    disp(rm.sourceData.expID{iExp})
    
    % Select and split up data from current experiment
    currModelData = rm.modelData(iExp, :);
%     currModelData.fullDataTbl{:} = currModelData.fullDataTbl{:}(:, [1 2 4]);
    tblFit = currModelData.fullDataTbl{:}(currModelData.fitRowInds{:}, :);
    tblTest = currModelData.fullDataTbl{:}(currModelData.testRowInds{:}, :);
    tblPred = currModelData.fullDataTbl{:}; 
    
    % Use all training data to create and evaluate a stepwise model
    kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
            mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
    emptyArgs = cellfun(@isempty, kvArgs);
    kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
    fullMdls{iExp} = stepwiselm(tblFit, kvArgs{:});
    
    fullMdls{iExp} = fitlm(tblFit, 'fl ~ 1 + fwSpeed + odorResp + fwSpeed:odorResp');

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


%% =================================================================================================
% PLOTTING RESULTS
% ==================================================================================================

%% Plot trial-averaged behavioral response to odor
saveFig = 0;
matchYLims = 1;
try
f = figure(2);clf;hold on
f.Color = [1 1 1];
f.Position = [f.Position(1:2), 1070, 950];
if (f.Position(2) + f.Position(4)) > 1000
    f.Position(2) = 1000 - f.Position(4);
end
allYlim = [];
for iFig = 1:size(rm.sourceData, 1)
    ax = subaxis(3, 3, iFig, 'ml', 0.06, 'mr', 0.03, 'sv', 0.06, 'mb', 0.06); hold on;
    rm.plot_mean_moveSpeed(iFig, ax);
    allYlim(iFig, :) = ylim();
    ax.Title.String = ax.Title.String(1:10);
    if ismember(iFig, [1 4 7])
        ylabel('Move speed (mm/sec)')
    end
    if iFig < 7
        xlabel(ax, '')
    end
    ax.FontSize = 12;
    box off
    xlim([0 6]); 
end
suptitle('Trial-averaged moveSpeed responses to odor stim');
f.Children(2).Children.FontWeight = 'bold';
f.Children(2).Children.FontSize = 16;
if matchYLims
    for iAx = [1 3:size(f.Children)]
        ylim(f.Children(iAx), [min(allYlim(:, 1)), max(allYlim(:, 2))])
    end
end
f.UserData.sourceDataParams = rm.sourceDataParams;
f.UserData.modelParams = rm.modelParams;
saveFileName = ['behavioral_odor_responses'];
if saveFig
    save_figure(f, figDir, saveFileName);
end
catch ME; rethrow(ME); end

%% Plot overlaid odor response filter and trial-averaged behavioral responses to odor
saveFig = 0;
try
f = figure(3);clf;hold on
f.Color = [1 1 1];
f.Position = [f.Position(1:2), 1070, 950];
if (f.Position(2) + f.Position(4)) > 1000
    f.Position(2) = 1000 - f.Position(4);
end
allYlim = [];
for iFig = 1:size(rm.sourceData, 1)
    ax = subaxis(3, 3, iFig, 'ml', 0.04, 'mr', 0.03, 'sv', 0.06, 'mb', 0.06); hold on;
    rm.plot_mean_moveSpeed(iFig, ax);
    allYlim(iFig, :) = ylim();
    yyaxis right
    rm.plot_odor_filter(iFig, ax);
    ax.Title.String = ax.Title.String(1:10);
    if iFig < 7
        xlabel(ax, '')
    end
    ax.YAxis(2).Visible = 'off';
    box off
    ax.YTick = [];
    xlim([0 6]);
    
end
for iPlot = 1:numel(f.Children)
    yyaxis(f.Children(iPlot), 'left')
    f.Children(iPlot).YLim = [min(allYlim(:, 1)), max(allYlim(:, 2))];
end
suptitle('Trial-averaged responses to odor stim (blue = GCaMP, black = move speed)');
f.Children(2).Children.FontWeight = 'bold';
f.Children(2).Children.FontSize = 16;
f.UserData.sourceDataParams = rm.sourceDataParams;
f.UserData.modelParams = rm.modelParams;
saveFileName = ['behavioral_and_estimated_odor_response_filters_overlay'];
if saveFig
    save_figure(f, figDir, saveFileName);
end
catch ME; rethrow(ME); end

%% Plot summary of adjusted R2 scores
saveFig = 0;
try
plotArr = [rm.modelData.fullMdlAdjR2, ...
        rm.modelData.driftCorrectedMdlAdjR2]';
    
groupNames = {'Base model', 'Drift corrected model'};
% plotArr = [rm.modelData.driftCorrectedMdlAdjR2, allRm{1}.modelData.driftCorrectedMdlAdjR2, ...
%         rm.modelData.noOdorHistMdlAdjR2]';

f = figure(3);clf; hold on 
f.Color = [1 1 1];
f.Position = [f.Position(1:2), 500, 500];
if (f.Position(2) + f.Position(4)) > 1000
    f.Position(2) = 1000 - f.Position(4);
end
plot(plotArr, '-o', 'Color', rgb('black'), ...
        'markerfacecolor', rgb('green'), 'markerEdgeColor', rgb('green')); 
SEM = std_err(plotArr, 2);
errorbar(1:size(plotArr, 1), mean(plotArr, 2)', SEM, '-s', ...
        'color', 'k', 'markerFaceColor', 'k', 'linewidth', 3, 'markerSize', 10)
xlim([0.5, size(plotArr, 1) + 0.5])
ylim([-0.001, 1])
ax = gca();
ax.XTick = 1:size(plotArr, 1);
ax.XTickLabel = groupNames;
ylabel('Model adjusted R^2');
ax.FontSize = 14;
if strcmpi(rm.modelParams.upper, 'linear')
    titleStr = 'No interaction terms';
    fileNameStr = 'no';
else
    titleStr = 'Interaction terms allowed';
    fileNameStr = 'with';
end
title({['Adj. R^2 step thresh: ', num2str(rm.modelParams.pEnter)], titleStr});
f.UserData.sourceDataParams = rm.sourceDataParams;
f.UserData.modelParams = rm.modelParams;
if saveFig
    pEnterStr = num2str(rm.modelParams.pEnter);
    saveFileName = ['adj_R2_summary_stepThresh_', pEnterStr(3:end), ...
            '_', fileNameStr, '_interactions'];
    save_figure(f, figDir, saveFileName);
end
catch ME; rethrow(ME); end

%% Plot summary grid of coefficients used in the drift-corrected models

try
allDcCoeffNames = {};
for iExp = 1:size(rm.modelData, 1)
    allDcCoeffNames = [allDcCoeffNames, rm.modelData.driftCorrectedMdls{iExp}.CoefficientNames];    
end
uniqueCoeffs = unique(allDcCoeffNames);
uniqueCoeffs = uniqueCoeffs(2:end) % Drop intercept term
%%
saveFig = 0;
uniqueCoeffs = uniqueCoeffs([1 3 2]);
coeffArrDc = zeros(numel(uniqueCoeffs), size(rm.modelData, 1));
for iExp = 1:size(rm.modelData, 1) 
    for iCoeff = 1:numel(uniqueCoeffs)
        if ismember(uniqueCoeffs{iCoeff}, rm.modelData.driftCorrectedMdls{iExp}.CoefficientNames)
            coeffArrDc(iCoeff, iExp) = rm.modelData.driftCorrectedMdls{iExp}.Coefficients.Estimate(...
                    strcmp(rm.modelData.driftCorrectedMdls{iExp}.CoefficientNames, ...
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
    text(1.5, iRow, ['   Adj. R^2 = ', num2str(rm.modelData.driftCorrectedMdlAdjR2(iRow), 2)], ...
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


predictorVars = {'fwSpeed'}; % odorResp, fwSpeed, yawSpeed
legendLocs = {'sw', 'sw', 'nw', 'nw'};
% legendLocs = {'best', 'best', 'best', 'best'};
% xLimits = {[], [], [300 600], [1800 3200], [700 1900], [300 1000], [300 600], [100 1100], [1000 1600]}';
xLimits = repmat({[]}, 1, 10);
% xLimits = {[400 800], [200 600], [1 500], [1 500], [], [200 600], [500 800], [1 400]}'; 

try
legendStr = {'predicted fl', 'measured fl'};
for iExp = 1:size(rm.modelData, 1)
    
    f = figure(iExp + 100); clf;
    f.Color = [1 1 1];
    
    currModelData = rm.modelData(iExp, :);
    currSourceData = rm.sourceData(iExp, :);
       
    fullDataTbl = currModelData.fullDataTbl{:};
    fullDataTbl = fullDataTbl(~isnan(fullDataTbl.fl), :);
    
    legendStr = {'predicted fl (drift-corrected model)', ...
            'measured fl'};
    dcPredFl = currModelData.driftCorrectedMdlPredFl{:};
    measuredFl = currModelData.driftCorrectedMdlMeasuredFl{:};
    xx = currSourceData.volTimes{:}(1:numel(measuredFl))';
    xL = xLimits{iExp};
    if isempty(xL)
        xL = [0 xx(end)];
    end
    
    % Measured Fl data
    ax = subaxis(24, 1, 7:16, 'mt', 0.05, 'mb', 0.07, 'ml', 0.03, 'mr', 0.03, 'sv', 0.01); hold on;
    plot(xx, measuredFl, 'linewidth', 1, 'color', rgb('green'));
%     plot(xx, dcPredFl, 'linewidth', 1, 'color', rgb('orangered'));
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
%     title([currModelData.expID{:}, '  —  drift-corrected model (adj. R^2 = ', ...
%             num2str(currModelData.driftCorrectedMdlAdjR2, 2), ')'])
    box off;
    
    % Drift-corrected model
    ax = subaxis(24, 1, 17:24, 'mt', 0.05, 'mb', 0.07, 'ml', 0.03, 'mr', 0.03, 'sv', 0.01); hold on;
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
    box off;
    
    % Predictor variables
    ax = subaxis(24, 1, 1:6); hold on
    legendStr = {};
    if any(strcmpi('odorResp', predictorVars))
        odorResp = fullDataTbl.odorResp(~logical(sum(isnan(table2array(fullDataTbl)), 2)));
        plot(xx, odorResp ./ max(abs([min(odorResp), max(odorResp)])), 'color', rgb('red'), ...
            'linewidth', 1);
        legendStr{end + 1} = 'Odor response';
    end
    if any(strcmpi('fwSpeed', predictorVars))
        plot(xx, fullDataTbl.fwSpeed(~logical(sum(isnan(table2array(fullDataTbl)), 2))), '-', ...
                'color', rgb('blue'), 'linewidth', 1);
        legendStr{end + 1} = 'abs(fwSpeed)';
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
            num2str(currModelData.driftCorrectedMdlAdjR2, 2), ')'])
        
    if saveFigs
        figTitle = ['predictedVsMeasuredFl_', fileNameStr, currModelData.expID{:}];
        save_figure(f, figDir, figTitle);
    end
end

catch ME; rethrow(ME); end

% ==================================================================================================
%% Show scatter plots of predictor variables vs. fluorescence
try
f = figure(102); clf;
f.Color = [1 1 1];
varNames = finalMdl.PredictorNames;
nPlots = numSubplots(numel(varNames));
if ~isempty(odorIntegrationWin)
    tblPred = tbl(:, [1:(odorHistoryCols(1) - 1), odorHistoryCols(bestWinSize), size(tbl, 2)]);
else
    tblPred = tbl;
end
tblPred = tblPred(~isnan(tblPred.fl), :);
suptitle('Predictor variables vs. measured fluorescence')
f.Children(end).Children.FontSize = 16;
f.Children(end).Children.FontWeight = 'bold';
for iPlot = 1:numel(varNames)
    ax = subaxis(nPlots(1), nPlots(2), iPlot, ...
            'mb', 0.08, 'ml', 0.07, 'mr', 0.04, 'sv', 0.08, 'sh', 0.07); 
    hold on;
    plot(tblPred.(varNames{iPlot}), tblPred.fl, '.', 'color', 'k')
    [p, S] = polyfit(tblPred.(varNames{iPlot}), tblPred.fl, 1);
    fit = polyval(p, tblPred.(varNames{iPlot}), 'omitnan');
    plt = plot(tblPred.(varNames{iPlot}), fit, '-', 'color', 'r', 'linewidth', 2.5);
    R = corrcoef(tblPred.(varNames{iPlot}), tblPred.fl);
    Rsquared = R(1,2)^2;
    legend(plt, [' y = ', num2str(round(p(1), 2)), 'x,  R^2 = ', num2str(Rsquared, 2)], ...
            'fontSize', 12)
    xlabel(regexprep(varNames{iPlot}, '_', '\\_'));
    ylabel('Fluorescence');
    ax.FontSize = 14;
end
catch ME; rethrow(ME); end
