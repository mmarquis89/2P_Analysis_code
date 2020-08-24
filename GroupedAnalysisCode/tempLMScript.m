parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
dataDir = fullfile(parentDir, 'all_experiments');

%% Load one or more existing model objects
threshStr = {'01', '015', '02'};
allRm = {};
for i=1:numel(threshStr)
    load(fullfile(parentDir, ['rm_stepThresh_', threshStr{i}, '_with_interactions.mat']), 'rm');
    allRm{i} = rm;
end

%% Set up a new model 

% Set source data parameters;
p = [];
p.roiName = 'TypeD';
p.maxSpeed = 100;
p.smWinVols = 5;
p.smWinFrames = 3;
p.smReps = 10;
p.ftLagVols = 3;
p.odorRespFilterDur = [7 5];

% Set model parameters
mp = [];
mp.trainTestSplit = 0.8;
mp.kFold = 100;
mp.criterion = 'adjrsquared'; % 'sse, 'aic', 'bic', '', or 'adjrsquared'
mp.upper = ['linear'];
mp.pEnter = [0.02];
mp.pRemove = [0];
mp.verbose = 0;
mp.odorIntegrationWin = [20:10:200];

expIDList = {'20190304-1', '20190315-1', '20190315-2', '20190315-3', '20190401-1', ... 
               '20190401-2', '20190411-1', '20190411-2', '20190411-3'};       

skipTrials = {[], [], [], [], [], ...
              [], [], [], []};
          
skipVols = {[], [1:1500], [], [], [], ...
              [], [], [], [1:2200]};
                     
try          
expInfoTbl = table(expIDList', skipTrials', skipVols', 'VariableNames', {'expID', 'skipTrials', ...
        'skipVols'});

% Create analysis object
rm = RegressionModelAnalysis(expInfoTbl, p);

% Initialize models
rm = rm.initialize_models(mp);

% Find optimal odor integration windows for each experiment
rm = rm.optimize_integration_windows();

catch ME; rethrow(ME); end

%% Use the calculated best window sizes to fit new models for each of the experiments 
try
fullMdls = {};
fullMdlPredFl = {};
fullMdlAdjR2 = [];
disp('Training final models...')
for iExp = 1:size(rm.sourceData, 1)
    if numel(rm.modelParams) > 1
        mp = rm.modelParams(iExp);
    else
        mp = rm.modelParams;
    end
    disp(rm.sourceData.expID{iExp})
    
    % Select data with best odorIntegration win for current experiment
    currModelData = rm.modelData(iExp, :);
    bestWinSize = currModelData.bestWinSizes;
    varNames = currModelData.fullDataTbl{:}.Properties.VariableNames;
    odorHistVarInds = ~cellfun(@isempty,regexp(varNames, 'odorHistory')); 
    bestHistVarInd = ~cellfun(@isempty, regexp(varNames, ['_', num2str(bestWinSize), '$']));
    tblFit = currModelData.fullDataTbl{:}(currModelData.fitRowInds{:}, ...
            logical(~odorHistVarInds + bestHistVarInd));
    tblTest = currModelData.fullDataTbl{:}(currModelData.testRowInds{:}, logical(~odorHistVarInds ...
            + bestHistVarInd));
    tblPred = currModelData.fullDataTbl{:}(:, logical(~odorHistVarInds + bestHistVarInd));

    % Use all training data to create and evaluate a final stepwise model
    kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
            mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
    emptyArgs = cellfun(@isempty, kvArgs);
    kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
    fullMdls{iExp} = stepwiselm(tblFit, kvArgs{:});
    [predFl, ~] = predict(fullMdls{iExp}, tblTest(~isnan(tblTest.fl), 1:end-1));
    [~, fullMdlAdjR2(iExp)] = r_squared(tblTest.fl(~isnan(tblTest.fl)), predFl, ...
            fullMdls{iExp}.NumCoefficients);
    [fullMdlPredFl{iExp}, ~] = predict(fullMdls{iExp}, tblPred(~isnan(tblPred.fl), 1:end-1));
end
disp('Final models created');
rm.modelData.fullMdls = fullMdls';
rm.modelData.fullMdlPredFl = fullMdlPredFl';
rm.modelData.fullMdlAdjR2 = fullMdlAdjR2';
disp(rm.modelData)

catch ME; rethrow(ME); end

%% Subtract volsFromOdorStart from fl data and re-fit models to "drift-corrected" data
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
    
    % Select data with best odorIntegration win for current experiment
    currModelData = rm.modelData(iExp, :);
    bestWinSize = currModelData.bestWinSizes;
    varNames = currModelData.fullDataTbl{:}.Properties.VariableNames;
    odorHistVarInds = ~cellfun(@isempty,regexp(varNames, 'odorHistory')); 
    bestHistVarInd = ~cellfun(@isempty, regexp(varNames, ['_', num2str(bestWinSize), '$']));
    tblPred = currModelData.fullDataTbl{:}(:, logical(~odorHistVarInds + bestHistVarInd));
        
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
    driftCorrectedMdlMeasuredFl{iExp} = tblPred.fl(~isnan(tblPred.fl));
    driftCorrectedMdlTestFl{iExp} = tblTest.fl(~isnan(tblTest.fl));
    
    % Use all training data to create and evaluate a final stepwise model
    kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
        mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
    emptyArgs = cellfun(@isempty, kvArgs);
    kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
    driftCorrectedMdls{iExp} = stepwiselm(tblFit, kvArgs{:});
    [predFl, ~] = predict(driftCorrectedMdls{iExp}, ...
            tblTest(~isnan(tblTest.fl), 1:end-1));
    [~, driftCorrectedMdlAdjR2(iExp)] = r_squared(driftCorrectedMdlTestFl{iExp}, ...
            predFl, driftCorrectedMdls{iExp}.NumCoefficients);
    [driftCorrectedMdlPredFl{iExp}, ~] = predict(driftCorrectedMdls{iExp}, ...
            tblPred(~isnan(tblPred.fl), 1:end-1));
end
disp('Final models created');
rm.modelData.driftCorrectedMdls = driftCorrectedMdls';
rm.modelData.driftCorrectedMdlMeasuredFl = driftCorrectedMdlMeasuredFl'; 
rm.modelData.driftCorrectedMdlTestFl = driftCorrectedMdlTestFl';
rm.modelData.driftCorrectedMdlPredFl = driftCorrectedMdlPredFl';
rm.modelData.driftCorrectedMdlAdjR2 = driftCorrectedMdlAdjR2';
disp(rm.modelData)

catch ME; rethrow(ME); end

%% Re-fit the drift-corrected models without the odor history variables
try
noOdorHistReFitMdls = {};
noOdorHistReFitMdlPredFl = {};
noOdorHistReFitMdlAdjR2 = [];
disp('Training no-odor history models...')
for iExp = 1:size(rm.sourceData, 1)
    if numel(rm.modelParams) > 1
        mp = rm.modelParams(iExp);
    else
        mp = rm.modelParams;
    end
    disp(rm.sourceData.expID{iExp})
    
    % Select data with best odorIntegration win for current experiment
    currModelData = rm.modelData(iExp, :);
    varNames = currModelData.fullDataTbl{:}.Properties.VariableNames;
    odorHistVarInds = ~cellfun(@isempty,regexp(varNames, 'odorHistory')); 
    tblPred = currModelData.fullDataTbl{:}(:, ~odorHistVarInds);
    
    % Subtract volsFromExpStart offset from fluorescence data
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
    
    % Use all training data to create and evaluate a final stepwise model
    kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
        mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
    emptyArgs = cellfun(@isempty, kvArgs);
    kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
    noOdorHistReFitMdls{iExp} = stepwiselm(tblFit, kvArgs{:});
    [predFl, ~] = predict(noOdorHistReFitMdls{iExp}, tblTest(~isnan(tblTest.fl), 1:end-1));
    [~, noOdorHistReFitMdlAdjR2(iExp)] = r_squared(rm.modelData.driftCorrectedMdlTestFl{iExp}, predFl, ...
            noOdorHistReFitMdls{iExp}.NumCoefficients);
    [noOdorHistReFitMdlPredFl{iExp}, ~] = predict(noOdorHistReFitMdls{iExp}, ...
           tblPred(~isnan(tblPred.fl), 1:end-1));
end
disp('Final models created');
rm.modelData.noOdorHistReFitMdls = noOdorHistReFitMdls';
rm.modelData.noOdorHistReFitMdlPredFl = noOdorHistReFitMdlPredFl';
rm.modelData.noOdorHistReFitMdlAdjR2 = noOdorHistReFitMdlAdjR2';
disp(rm.modelData)
catch ME; rethrow(ME); end

%% Remove the odor history variables from the drift-corrected models
try
noOdorHistMdls = {};
noOdorHistMdlPredFl = {};
noOdorHistMdlAdjR2 = [];
disp('Training no-odor history models...')
for iExp = 1:size(rm.sourceData, 1)
    if numel(rm.modelParams) > 1
        mp = rm.modelParams(iExp);
    else
        mp = rm.modelParams;
    end
    disp(rm.sourceData.expID{iExp})
    
    % Select data with best odorIntegration win for current experiment
    currModelData = rm.modelData(iExp, :);
    varNames = currModelData.fullDataTbl{:}.Properties.VariableNames;
    odorHistVarInds = ~cellfun(@isempty,regexp(varNames, 'odorHistory')); 
    tblPred = currModelData.fullDataTbl{:}(:, ~odorHistVarInds);
    
    % Subtract volsFromExpStart offset from fluorescence data
    currMdl = rm.modelData.fullMdls{iExp};
    if ismember('volsFromExpStart', currMdl.Coefficients.Properties.RowNames)
        tblPredVolAdjust = tblPred;
        tblPredVolAdjust.fl = tblPredVolAdjust.fl - (currMdl.Coefficients{'volsFromExpStart', ...
                'Estimate'} .* tblPredVolAdjust.volsFromExpStart);
        tblPred = tblPredVolAdjust(:, ~strcmp(tblPredVolAdjust.Properties.VariableNames, ...
                'volsFromExpStart'));
    end
    tblTest = tblPred(currModelData.testRowInds{:}, :);
    
    % Remove all terms involving odor history from the model
    currModelData = rm.modelData(iExp, :);
    driftCorrectedMdl = currModelData.driftCorrectedMdls{:};
    coeffNames = driftCorrectedMdl.CoefficientNames;
    odorHistCoeffs = coeffNames(~cellfun(@isempty, regexp(coeffNames, 'odorHistory')));
    updatedMdl = driftCorrectedMdl;
    for iCoeff = 1:numel(odorHistCoeffs)
            updatedMdl = removeTerms(updatedMdl, odorHistCoeffs{iCoeff});
    end
    noOdorHistMdls{iExp} = updatedMdl;
    
    % Evaluate updated model
    [predFl, ~] = predict(noOdorHistMdls{iExp}, tblTest(~isnan(tblTest.fl), ...
            ~ismember(tblTest.Properties.VariableNames, odorHistCoeffs)));
    [~, noOdorHistMdlAdjR2(iExp)] = r_squared(tblTest.fl(~isnan(tblTest.fl)), predFl, ...
            noOdorHistMdls{iExp}.NumCoefficients);
    [noOdorHistMdlPredFl{iExp}, ~] = predict(noOdorHistMdls{iExp}, ...
           tblPred(~isnan(tblPred.fl), 1:end-1));
end
disp('Final models created');
rm.modelData.noOdorHistMdls = noOdorHistMdls';
rm.modelData.noOdorHistMdlPredFl = noOdorHistMdlPredFl';
rm.modelData.noOdorHistMdlAdjR2 = noOdorHistMdlAdjR2';
disp(rm.modelData)

catch ME; rethrow(ME); end

%% Plot summary of adjusted R2 scores
try
plotArr = [rm.modelData.fullMdlAdjR2, rm.modelData.driftCorrectedMdlAdjR2, ...
        rm.modelData.noOdorHistMdlAdjR2]';
f = figure(6);clf; hold on 
f.Color = [1 1 1];
plot(plotArr, '-o', 'Color', rgb('black'), ...
        'markerfacecolor', rgb('green'), 'markerEdgeColor', rgb('green')); 
SEM = std_err(plotArr, 2);
errorbar(1:size(plotArr, 1), mean(plotArr, 2)', SEM, '-s', ...
        'color', 'k', 'markerFaceColor', 'k', 'linewidth', 3, 'markerSize', 10)
xlim([0.5, size(plotArr, 1) + 0.5])
ylim([0, 1])
ax = gca();
ax.XTick = 1:size(plotArr, 1);
ax.XTickLabel = {'Initial model', 'Drift-corrected', 'No odor history', 'No odor history (re-fit)'};
ylabel('Model adjusted R^2');
ax.FontSize = 14;
title({['Adj. R^2 step thresh: ', num2str(rm.modelParams.pEnter)], 'Interaction terms allowed'});
catch ME; rethrow(ME); end

%% Plot all model coefficients
try
    
% Full model
f = figure(1); clf; hold on;
f.Color = [1 1 1];
subplotDims = [2 5];
for iExp = 1:size(rm.modelData, 1)
    
    currMdl = rm.modelData.fullMdls{iExp};
    
    subaxis(subplotDims(1), subplotDims(2), iExp, 'mt', 0.05, 'mb', 0.12, 'sv', 0.16, ...
            'ml', 0.05, 'mr', 0.05);
    hold on
    ax = gca();
    RegressionModelAnalysis.plot_coeffs(currMdl, ax);
    title([regexprep(rm.modelData.expID{iExp}, '_', '\\_'), ' coeffs (adj. R^2 = ', ...
            num2str(rm.modelData.fullMdlAdjR2(iExp), 2), ')'])
end

% Drift-corrected model
f = figure(2); clf; hold on;
f.Color = [1 1 1];
subplotDims = [2 5];
for iExp = 1:size(rm.modelData, 1)
    currMdl = rm.modelData.driftCorrectedMdls{iExp};
    subaxis(subplotDims(1), subplotDims(2), iExp, 'mt', 0.05, 'mb', 0.12, 'sv', 0.16, ...
            'ml', 0.05, 'mr', 0.05);
    hold on
    ax = gca();
    RegressionModelAnalysis.plot_coeffs(currMdl, ax);
    title([regexprep(rm.modelData.expID{iExp}, '_', '\\_'), ' coeffs (adj. R^2 = ', ...
            num2str(rm.modelData.driftCorrectedMdlAdjR2(iExp), 2), ')'])
end

% No-odor history model
f = figure(3); clf; hold on;
f.Color = [1 1 1];
subplotDims = [2 5];
for iExp = 1:numel(noOdorHistMdls)
    currMdl = noOdorHistMdls{iExp};
    if currMdl.NumCoefficients > 1
        subaxis(subplotDims(1), subplotDims(2), iExp, 'mt', 0.05, 'mb', 0.12, 'sv', 0.16, ...
            'ml', 0.05, 'mr', 0.05);
        hold on
        ax = gca();
        RegressionModelAnalysis.plot_coeffs(currMdl, ax);
        title([regexprep(rm.modelData.expID{iExp}, '_', '\\_'), ' coeffs (adj. R^2 = ', ...
            num2str(noOdorHistMdlAdjR2(iExp), 2), ')'])
    end
end
catch ME; rethrow(ME); end

%% Plot predicted vs. measured fluorescence
try
predictorVars = {'odorHistory'}; % odorHistory, odorResp, fwSpeed, yawSpeed

legendStr = {'predicted fl', 'measured fl'};
for iExp = 1:size(rm.modelData, 1)
    
    f = figure(iExp + 10); clf;
    f.Color = [1 1 1];
    
    currModelData = rm.modelData(iExp, :);
    currSourceData = rm.sourceData(iExp, :);
    
    bestWinSize = currModelData.bestWinSizes;
    
    fullDataTbl = currModelData.fullDataTbl{:};
    fullDataTbl = fullDataTbl(~isnan(fullDataTbl.fl), :);
    
    
    % Full model
    ax = subaxis(2, 1, 1, 'mt', 0.05, 'mb', 0.07, 'ml', 0.03, 'mr', 0.03); 
    hold on
%     legendStr = {'predicted fl', 'measured fl'};
%     predFl = currModelData.fullMdlPredFl{:}; 
%     xx = currSourceData.volTimes{:}(1:numel(predFl))';
%     plot(xx, predFl, 'linewidth', 1, 'color', rgb('orangered'));
%     plot(xx, fullDataTbl.fl, 'linewidth', 1, 'color', rgb('darkmagenta'));
%     ax.YTickLabel = [];
%     yyaxis right; hold on
%     if any(strcmpi('odorResp', predictorVars))
%         plot(xx, fullDataTbl.odorResp, 'color', rgb('blue'));
%         legendStr{end + 1} = 'Odor response';
%     end
%     if any(strcmpi('fwSpeed', predictorVars))
%         plot(xx, fullDataTbl.fwSpeed, '-', 'color', rgb('gold'));
%         legendStr{end + 1} = 'abs(fwSpeed)';
%     end
%     if any(strcmpi('yawSpeed', predictorVars))
%         plot(xx, fullDataTbl.yawSpeed, 'color', rgb('red'));
%         legendStr{end + 1} = 'yawSpeed';
%     end
%     if any(strcmpi('odorHistory', predictorVars))
%         varName = ['odorHistory_', num2str(bestWinSize)];
%         plot(xx, fullDataTbl.(varName), '-', 'color', rgb('green'));
%         legendStr{end + 1} = regexprep(varName, '_', '\\_');
%     end
%     legend(legendStr, 'fontsize', 12, 'Location', 'best', 'autoupdate', 'off')
%     ax.FontSize = 14;
%     ax.YTickLabel = [];
%     ax.XTickLabel = [];
%     ax.XLim = [0, xx(end)];
%     title([currModelData.expID{:}, '  —  Full model (adj. R^2 = ', ...
%             num2str(currModelData.fullMdlAdjR2, 2), ')'])
    
    % Drift corrected model
%     ax = subaxis(3, 1, 2); hold on
    legendStr = {'predicted fl', 'measured fl'};
    predFl = currModelData.driftCorrectedMdlPredFl{:};
    measuredFl = currModelData.driftCorrectedMdlMeasuredFl{:};
    xx = currSourceData.volTimes{:}(1:numel(predFl))';
    plot(xx, predFl, 'linewidth', 1, 'color', rgb('orangered'));
    plot(xx, measuredFl, 'linewidth', 1, 'color', rgb('darkmagenta'));
    ax.YTickLabel = [];
    yyaxis right; hold on
    if any(strcmpi('odorResp', predictorVars))
        plot(xx, fullDataTbl.odorResp, 'color', rgb('blue'));
        legendStr{end + 1} = 'Odor response';
    end
    if any(strcmpi('fwSpeed', predictorVars))
        plot(xx, fullDataTbl.fwSpeed, '-', 'color', rgb('gold'));
        legendStr{end + 1} = 'abs(fwSpeed)';
    end
    if any(strcmpi('yawSpeed', predictorVars))
        plot(xx, fullDataTbl.yawSpeed, 'color', rgb('red'));
        legendStr{end + 1} = 'yawSpeed';
    end
    if any(strcmpi('odorHistory', predictorVars))
        varName = ['odorHistory_', num2str(bestWinSize)];
        plot(xx, fullDataTbl.(varName), '-', 'color', rgb('green'));
        legendStr{end + 1} = regexprep(varName, '_', '\\_');
    end
    legend(legendStr, 'fontsize', 12, 'Location', 'best', 'autoupdate', 'off');
    ax.FontSize = 14;
    ax.YTickLabel = [];
    ax.XTickLabel = [];
    ax.XLim = [0, xx(end)];
    title([currModelData.expID{:}, '  —  drift-corrected model (adj. R^2 = ', ...
            num2str(currModelData.driftCorrectedMdlAdjR2, 2), ')'])
    
    % No-odor history model
    ax = subaxis(2, 1, 2); hold on
    legendStr = {'predicted fl', 'measured fl'};
    predFl = currModelData.noOdorHistMdlPredFl{:};
    measuredFl = currModelData.driftCorrectedMdlMeasuredFl{:};
    xx = currSourceData.volTimes{:}(1:numel(predFl))';
    plot(xx, predFl, 'linewidth', 1, 'color', rgb('orangered'));
    plot(xx, measuredFl, 'linewidth', 1, 'color', rgb('darkmagenta'));
    ax.YTickLabel = [];
    yyaxis right; hold on
    if any(strcmpi('odorResp', predictorVars))
        plot(xx, fullDataTbl.odorResp, 'color', rgb('blue'));
        legendStr{end + 1} = 'Odor response';
    end
    if any(strcmpi('fwSpeed', predictorVars))
        plot(xx, fullDataTbl.fwSpeed, '-', 'color', rgb('gold'));
        legendStr{end + 1} = 'abs(fwSpeed)';
    end
    if any(strcmpi('yawSpeed', predictorVars))
        plot(xx, fullDataTbl.yawSpeed, 'color', rgb('red'));
        legendStr{end + 1} = 'yawSpeed';
    end
    if any(strcmpi('odorHistory', predictorVars))
        varName = ['odorHistory_', num2str(bestWinSize)];
        plot(xx, fullDataTbl.(varName), '-', 'color', rgb('green'));
        legendStr{end + 1} = regexprep(varName, '_', '\\_');
    end
    legend(legendStr, 'fontsize', 12, 'Location', 'best', 'autoupdate', 'off');
    ax.FontSize = 14;
    xlabel('Time (sec)')
    ax.YTickLabel = [];
    ax.XLim = [0, xx(end)];
    title(['Odor history-free model (adj. R^2 = ', ...
            num2str(currModelData.noOdorHistMdlAdjR2, 2), ')'])
    
end

catch ME; rethrow(ME); end




%% Combine data 
try
allExpDataTbl = [];
for iExp = 1:size(rm.modelData, 1)
    
    allExpDataTbl = [allExpDataTbl; rm.modelData.fullDataTbl{iExp}];
    
end
shuffleInds = randperm(size(allExpDataTbl, 1));
splitInd = floor(numel(shuffleInds) * mp.trainTestSplit);
fitRowInds = shuffleInds(1:splitInd);
testRowInds = shuffleInds((splitInd + 1):end);


disp('Finding best odor integration window sizes...')

% Get training data
tblTrain = allExpDataTbl(fitRowInds, :);
    
% Drop rows without valid fluorescence measurements
tblTrain = tblTrain(~isnan(tblTrain.fl), :);
   
% Divide training data into cross-validation folds
chunkSize = floor(size(tblTrain, 1) / mp.kFold);
startInds = 1:chunkSize:size(tblTrain, 1);
endInds = [startInds(2:end) - 1, size(tblTrain, 1)];

% Identify positions of odor history variables
kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
    mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
emptyArgs = cellfun(@isempty, kvArgs);
kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
varNames = tblTrain.Properties.VariableNames;
odorHistVars = ~cellfun(@isempty,regexp(varNames, 'odorHistory'));
odorHistVarInds = find(odorHistVars);
    
% Select odor history window size using cross-validation
nWindows = numel(mp.odorIntegrationWin);
loopTbls = [];
for iWin = 1:nWindows
    loopTbls{iWin} = tblTrain(:, [1:(odorHistVarInds(1) - 1), odorHistVarInds(iWin), ...
        size(tblTrain, 2)]);
end
kFold = mp.kFold;
predAdjR2 = zeros(nWindows, kFold);
allCvMdls = cell(nWindows, kFold);
parfor iWin = 1:nWindows

    currTbl = loopTbls{iWin};

    % Fit model with cross-validation
    for iFold = 1:kFold
        if iWin == 1 && ~mod(iFold, 20)
            disp(['Running cross-validation fold ', num2str(iFold), ' of ', ...
                num2str(kFold)]);
        end

        % Split data for current fold
        currFoldTrainTbl = currTbl;
        currFoldTrainTbl(startInds(iFold):endInds(iFold), :) = [];
        currFoldTestTbl = currTbl(startInds(iFold):endInds(iFold), :);

        % ------------- Generate model ---------------
        mdl = stepwiselm(currFoldTrainTbl, kvArgs{:});
        % mdl = fitlm(tblFit);
        allCvMdls{iWin, iFold} = compact(mdl);

        % Evaluate fit to test data
        [noOdorHistMdlPredFl, ~] = predict(mdl, currFoldTestTbl(:, 1:end-1));
        [~, R2] = r_squared(currFoldTestTbl.fl, noOdorHistMdlPredFl, mdl.NumCoefficients);
        predAdjR2(iWin, iFold) = R2;

    end%iFold
end%iWin

bestWinSize = mp.odorIntegrationWin(argmax(mean(predAdjR2, 2), 1));
disp('Optimal odor integration window selected');

% Use the calculated best window sizes to fit new models for each of the experiments 

disp('Training final model...')

    
% Select data with best odorIntegration win for current experiment
odorHistVarInds = ~cellfun(@isempty,regexp(varNames, 'odorHistory')); 
bestHistVarInd = ~cellfun(@isempty, regexp(varNames, ['_', num2str(bestWinSize)]));
tblFit = allExpDataTbl(fitRowInds, logical(~odorHistVarInds + bestHistVarInd));
tblTest = allExpDataTbl(testRowInds, logical(~odorHistVarInds + bestHistVarInd));
    
% Use all training data to create and evaluate a final stepwise model
kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
        mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
emptyArgs = cellfun(@isempty, kvArgs);
kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
finalMdl= stepwiselm(tblFit, kvArgs{:});
[noOdorHistMdlPredFl, ~] = predict(finalMdl, tblTest(~isnan(tblTest.fl), 1:end-1));
[~, predAdjR2] = r_squared(tblTest.fl(~isnan(tblTest.fl)), noOdorHistMdlPredFl, ...
        finalMdl.NumCoefficients);

disp('Final model created');


catch ME; rethrow(ME); end

%% Show scatter plots of predictor variables vs. fluorescence
try
hdls{4} = figure(102); clf;
hdls{4}.Color = [1 1 1];
varNames = finalMdl.PredictorNames;
nPlots = numSubplots(numel(varNames));
if ~isempty(odorIntegrationWin)
    tblPred = tbl(:, [1:(odorHistoryCols(1) - 1), odorHistoryCols(bestWinSize), size(tbl, 2)]);
else
    tblPred = tbl;
end
tblPred = tblPred(~isnan(tblPred.fl), :);
suptitle('Predictor variables vs. measured fluorescence')
hdls{4}.Children(end).Children.FontSize = 16;
hdls{4}.Children(end).Children.FontWeight = 'bold';
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
