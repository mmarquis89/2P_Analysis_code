%% Set up a new model object

% Set source data parameters;
p = [];
p.maxSpeed = 100;
p.smWinVols = 5;
p.roiName = 'EB-DAN';
p.smWinFrames = 7;
p.smReps = 20;
p.ftLagVols = 3;
p.speedType = 'forwardSpeed';

expIDList = {'20201015-1', '20201015-2', '20201019-1', '20201019-2', '20201023-1', '20201023-2', ...
        '20201027-1', '20201029-1', '20201029-2', '20201029-3', '20201029-4'};

skipTrials = {[11:15], [11:16], [8:15], [1:4, 9, 12:13], 7:10, 8:12, 6:8, 7:10, 8:10, 7:9, 7:11};

skipVols = repmat({[]}, 1, numel(skipTrials));
                     
         
expInfoTbl = table(expIDList', skipTrials', skipVols', 'VariableNames', {'expID', 'skipTrials', ...
        'skipVols'});

% Create analysis object
rm = RegressionModelAnalysis_PPM3(expInfoTbl, p);

%% Choose parameters and initialize model
mp = [];
mp.trainTestSplit = 0.8;
mp.criterion = 'adjrsquared'; % 'sse, 'aic', 'bic', '', or 'adjrsquared'
mp.upper = [];
mp.pEnter = [0.03];
mp.pRemove = [0];
mp.verbose = 0;
mp.useYaw = 1;
mp.useRunSpeed = 1;
mp.useDriftCorrection = 0;

% Initialize models
rm = rm.initialize_models(mp);

%% Train initial models with and without yaw and running speed
try

mp.useYaw = 1;
mp.useRunSpeed = 1;
rm_full = rm.initialize_models(mp);
rm_full = rm_full.train_initial_models();

mp.useYaw = 0;
mp.useRunSpeed = 1;
rm_noYawSpeed = rm.initialize_models(mp);
rm_noYawSpeed.modelData.fitRowInds = rm_full.modelData.fitRowInds;
rm_noYawSpeed.modelData.testRowInds = rm_full.modelData.testRowInds;
rm_noYawSpeed = rm_noYawSpeed.train_initial_models(); 

mp.useYaw = 1;
mp.useRunSpeed = 0;
rm_noRunSpeed = rm.initialize_models(mp);
rm_noRunSpeed.modelData.fitRowInds = rm_full.modelData.fitRowInds;
rm_noRunSpeed.modelData.testRowInds = rm_full.modelData.testRowInds;
rm_noRunSpeed = rm_noRunSpeed.train_initial_models();    

rm = rm_full;

catch ME; rethrow(ME); end


%% Plot summary of adjusted R2 values for each condition

saveFig = 0;
try
    
plotArr = [rm_full.modelData.fullMdlAdjR2, rm_noRunSpeed.modelData.fullMdlAdjR2, ...
        rm_noYawSpeed.modelData.fullMdlAdjR2]';

% groupNames = {'Full model', ['No ', rm_noRunSpeed.sourceDataParams.speedType], 'No yawSpeed'};
groupNames = {'Full model', 'No fwSpeed', 'No yawSpeed'};

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
if strcmpi(rm_full.modelParams.upper, 'linear')
    titleStr = 'No interaction terms';
    fileNameStr = 'no';
else
    titleStr = 'Interaction terms allowed';
    fileNameStr = 'with';
end
title({['Adj. R^2 step thresh: ', num2str(rm_full.modelParams.pEnter)], titleStr});
f.UserData.sourceDataParams = rm_full.sourceDataParams;
f.UserData.modelParams = rm_full.modelParams;
if saveFig
    pEnterStr = num2str(rm_full.modelParams.pEnter);
    saveFileName = [p.roiName, '_adj_R2_summary_stepThresh_', pEnterStr(3:end), ...
            '_', fileNameStr, '_interactions'];
    save_figure(f, figDir, saveFileName);
end
    
    
    
        
catch ME; rethrow(ME); end


%% Plot summary grid of coefficients used in the models

saveFig = 0;
try
    
allCoeffNames = {};
for iExp = 1:size(rm.modelData, 1)
    allCoeffNames = [allCoeffNames, rm.modelData.fullMdls{iExp}.CoefficientNames];    
end
uniqueCoeffs = unique(allCoeffNames);
uniqueCoeffs = uniqueCoeffs(2:end); % Drop intercept term
%

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
    text(0.5, iRow, ['   Adj. R^2 = ', num2str(rm.modelData.fullMdlAdjR2(iRow), 2)], ...
            'FontSize', 12);
end
title('Model coefficients')
f.UserData.sourceDataParams = rm.sourceDataParams;
f.UserData.modelParams = rm.modelParams;
if saveFig
    pEnterStr = num2str(rm.modelParams.pEnter);
    saveFileName = ['FB-DAN_model_coeff_summary_fwSpeed_noInteractions'];
    save_figure(f, figDir, saveFileName);
end
catch ME; rethrow(ME); end


%% Plot predicted vs. measured fluorescence

saveFigs = 0;
fileNameStr = 'fullExp_';

predictorVars = {'fwSpeed', 'yawSpeed'}; % odorResp, fwSpeed, moveSpeed, yawSpeed
legendLocs = {'sw', 'sw', 'nw', 'nw'};
% legendLocs = {'best', 'best', 'best', 'best'};
xLimits = repmat({[]}, 1, size(rm.modelData,1));

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
        plot(xx, fullDataTbl.forwardSpeed(~logical(sum(isnan(table2array(fullDataTbl)), 2))), '-', ...
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