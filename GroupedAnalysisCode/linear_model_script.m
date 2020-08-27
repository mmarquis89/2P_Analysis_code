parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
dataDir = fullfile(parentDir, 'all_experiments');
figDir = ...
        'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Figs\linear_regression_analysis';

%% Load one or more existing model objects
threshStr = {'03'};%{'01', '015', '02'};
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
mp.upper = [];
mp.pEnter = [0.01];
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
%% =================================================================================================
%  MODEL POSTPROCESSING
%  =================================================================================================

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

%% =================================================================================================
% PLOTTING RESULTS
% ==================================================================================================

%% Plot estimated odor response filters
saveFig = 0;
try
f = figure(1);clf;hold on
f.Color = [1 1 1];
f.Position = [f.Position(1:2), 1070, 950];
if (f.Position(2) + f.Position(4)) > 1000
    f.Position(2) = 1000 - f.Position(4);
end
for iFig = 1:size(rm.sourceData, 1)
    ax = subaxis(3, 3, iFig, 'ml', 0.04, 'mr', 0.03, 'sv', 0.06, 'mb', 0.06); hold on;
    rm.plot_odor_filter(iFig, ax);
%     rm.plot_mean_moveSpeed(iFig, ax);
    ax.Title.String = ax.Title.String(1:10);
    if iFig < 7
        xlabel(ax, '')
    end
    box off
    ax.YTick = 0;
    xL = xlim();
    yL = ylim();
    if yL(1) ~= 0
        plot(xL, [0 0], '--', 'color', 'blue')
    end
end
suptitle('Estimated odor response filters');
f.Children(2).Children.FontWeight = 'bold';
f.Children(2).Children.FontSize = 16;
saveFileName = ['estimated_odor_response_filters'];
if saveFig
    save_figure(f, figDir, saveFileName);
end
catch ME; rethrow(ME); end

%% Plot trial-averaged behavioral response to odor
saveFig = 0;
try
f = figure(2);clf;hold on
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
%     yyaxis right
%     rm.plot_odor_filter(iFig, ax)
    ax.Title.String = ax.Title.String(1:10);
    if iFig < 7
        xlabel(ax, '')
    end
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
saveFileName = ['behavioral_and_GCaMP_odor_responses_overlay'];
if saveFig
    save_figure(f, figDir, saveFileName);
end
catch ME; rethrow(ME); end

%% Plot summary of adjusted R2 scores
saveFig = 1;
try
plotArr = [rm.modelData.driftCorrectedMdlAdjR2, ...
        rm.modelData.noOdorHistMdlAdjR2]';
f = figure(3);clf; hold on 
f.Color = [1 1 1];
f.Position = [f.Position(1:2), 570, 500];
if (f.Position(2) + f.Position(4)) > 1000
    f.Position(2) = 1000 - f.Position(4);
end
plot(plotArr, '-o', 'Color', rgb('black'), ...
        'markerfacecolor', rgb('green'), 'markerEdgeColor', rgb('green')); 
SEM = std_err(plotArr, 2);
errorbar(1:size(plotArr, 1), mean(plotArr, 2)', SEM, '-s', ...
        'color', 'k', 'markerFaceColor', 'k', 'linewidth', 3, 'markerSize', 10)
xlim([0.5, size(plotArr, 1) + 0.5])
ylim([0, 1])
ax = gca();
ax.XTick = 1:size(plotArr, 1);
ax.XTickLabel = {'Drift-corrected', 'No odor history', 'No odor history (re-fit)'};
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
if saveFig
    pEnterStr = num2str(rm.modelParams.pEnter);
    saveFileName = ['adj_R2_summary_stepThresh_', pEnterStr(3:end), ...
            '_', fileNameStr, '_interactions'];
    save_figure(f, figDir, saveFileName);
end
catch ME; rethrow(ME); end

%% Plot all model coefficients
saveFig = 0;
try
if strcmpi(rm.modelParams.upper, 'linear')
    titleStr = 'no interaction terms';
    fileNameStr = 'no';
else
    titleStr = 'interaction terms allowed';
    fileNameStr = 'with';
end

% Full model
f_1 = figure(1); clf; hold on;
f_1.Color = [1 1 1];
f_1.Position = [f_1.Position(1:2), 1900, 950];
if (f_1.Position(2) + f_1.Position(4)) > 1000
    f_1.Position(2) = 1000 - f_1.Position(4);
end
subplotDims = [2 5];
for iExp = 1:size(rm.modelData, 1)
    
    currMdl = rm.modelData.fullMdls{iExp};
    
    subaxis(subplotDims(1), subplotDims(2), iExp, 'mt', 0.05, 'mb', 0.12, 'sv', 0.16, ...
            'ml', 0.05, 'mr', 0.05);
    hold on
    ax = gca();
    RegressionModelAnalysis.plot_coeffs(currMdl, ax);
    title([regexprep(rm.modelData.expID{iExp}, '_', '\\_'), ' (adj. R^2 = ', ...
            num2str(rm.modelData.fullMdlAdjR2(iExp), 2), ')'], 'fontweight', 'normal')
end
suptitle(['Full model coefficients (stepThresh: ', num2str(rm.modelParams.pEnter), ', ', ...
        titleStr, ')']);
f_1.Children(2).Children.FontWeight = 'bold';
f_1.Children(2).Children.FontSize = 16;

% Drift-corrected model
f_2 = figure(2); clf; hold on;
f_2.Color = [1 1 1];
f_2.Position = [f_2.Position(1:2), 1900, 950];
if (f_2.Position(2) + f_2.Position(4)) > 1000
    f_2.Position(2) = 1000 - f_2.Position(4);
end
subplotDims = [2 5];
for iExp = 1:size(rm.modelData, 1)
    currMdl = rm.modelData.driftCorrectedMdls{iExp};
    subaxis(subplotDims(1), subplotDims(2), iExp, 'mt', 0.05, 'mb', 0.12, 'sv', 0.16, ...
            'ml', 0.05, 'mr', 0.05);
    hold on
    ax = gca();
    RegressionModelAnalysis.plot_coeffs(currMdl, ax);
    title([regexprep(rm.modelData.expID{iExp}, '_', '\\_'), ' (adj. R^2 = ', ...
            num2str(rm.modelData.driftCorrectedMdlAdjR2(iExp), 2), ')'], 'fontweight', 'normal')
end
suptitle(['Drift-corrected model coefficients (stepThresh: ', num2str(rm.modelParams.pEnter), ...
        ', ', titleStr, ')']);
f_2.Children(2).Children.FontWeight = 'bold';
f_2.Children(2).Children.FontSize = 16;

% No-odor history model
f_3 = figure(3); clf; hold on;
f_3.Color = [1 1 1];
f_3.Position = [f_3.Position(1:2), 1900, 950];
if (f_3.Position(2) + f_3.Position(4)) > 1000
    f_3.Position(2) = 1000 - f_3.Position(4);
end
subplotDims = [2 5];
for iExp = 1:numel(rm.modelData.noOdorHistMdls)
    currMdl = rm.modelData.noOdorHistMdls{iExp};
    if currMdl.NumCoefficients > 1
        subaxis(subplotDims(1), subplotDims(2), iExp, 'mt', 0.05, 'mb', 0.12, 'sv', 0.16, ...
            'ml', 0.05, 'mr', 0.05);
        hold on
        ax = gca();
        RegressionModelAnalysis.plot_coeffs(currMdl, ax);
        title([regexprep(rm.modelData.expID{iExp}, '_', '\\_'), ' (adj. R^2 = ', ...
            num2str(noOdorHistMdlAdjR2(iExp), 2), ')'], 'fontweight', 'normal')
    end
end
suptitle(['No-odor history model coefficients (stepThresh: ', num2str(rm.modelParams.pEnter), ...
        ', ', titleStr, ')']);
f_3.Children(2).Children.FontWeight = 'bold';
f_3.Children(2).Children.FontSize = 16;
if saveFig
    pEnterStr = num2str(rm.modelParams.pEnter);
    
    saveFileName = ['full_model_coeffs_stepThresh_', pEnterStr(3:end), '_', fileNameStr, ...
            '_interactions'];
    save_figure(f_1, figDir, saveFileName);
    saveFileName = ['driftCorrected_model_coeffs_stepThresh_', pEnterStr(3:end), '_', fileNameStr, ...
            '_interactions'];
    save_figure(f_2, figDir, saveFileName);
    saveFileName = ['noOdorHist_model_coeffs_stepThresh_', pEnterStr(3:end), '_', fileNameStr, ...
            '_interactions'];
    save_figure(f_3, figDir, saveFileName);
end
catch ME; rethrow(ME); end

%% Plot summary grid of coefficients used in the drift-corrected models
saveFig = 0;
try
allDcCoeffNames = {};
for iExp = 1:size(rm.modelData, 1)
    allDcCoeffNames = [allDcCoeffNames, rm.modelData.driftCorrectedMdls{iExp}.CoefficientNames];    
end
uniqueCoeffs = unique(regexprep(allDcCoeffNames, 'odorHistory_.*', 'odorHistory'));
uniqueCoeffs = uniqueCoeffs(2:end); % Drop intercept term
uniqueCoeffs = uniqueCoeffs([1 4 3 2 5]);
coeffArrDc = zeros(numel(uniqueCoeffs), size(rm.modelData, 1));
for iExp = 1:size(rm.modelData, 1)
%     coeffArrDc(ismember(uniqueCoeffs, regexprep(...
%             rm.modelData.driftCorrectedMdls{iExp}.CoefficientNames, 'odorHistory_.*', ...
%             'odorHistory')), iExp) = 1;    
    for iCoeff = 1:numel(uniqueCoeffs)
        if ismember(uniqueCoeffs{iCoeff}, ...
                regexprep(rm.modelData.driftCorrectedMdls{iExp}.CoefficientNames, ...
                'odorHistory_.*', 'odorHistory'))
            coeffArrDc(iCoeff, iExp) = rm.modelData.driftCorrectedMdls{iExp}.Coefficients.Estimate(...
                    strcmp(regexprep(rm.modelData.driftCorrectedMdls{iExp}.CoefficientNames, ...
                    'odorHistory_.*', 'odorHistory'), uniqueCoeffs{iCoeff}));
        end
    end
end
% coeffTbl = array2table(coeffArr, 'rowNames', uniqueCoeffs, 'variableNames', rm.modelData.expID);
coeffTblDc = array2table(coeffArrDc', 'variableNames', uniqueCoeffs, 'rowNames', ...
        rm.modelData.expID);

f = figure(35); clf; 
f.Color = [1 1 1];
f.Position = [f.Position(1:2), 1035, 985];
if (f.Position(2) + f.Position(4)) > 1000
    f.Position(2) = 1000 - f.Position(4);
end
plotArr = table2array(coeffTblDc);
lim = max(abs([max(plotArr(:)), min(plotArr(:))]));
xSize = size(plotArr, 2);
plotArr(:, end + 1) = -lim;
plotArr(:, end + 1) = lim;
imagesc(plotArr);
cm = colormap(bluewhitered);
xlim([0.5, xSize + 0.5])
ax = gca;
% ax.XTickLabel = rm.modelData.expID;
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
colorbar();
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
    text(4.5, iRow, ['   Adj. R^2 = ', num2str(rm.modelData.driftCorrectedMdlAdjR2(iRow), 2)], ...
            'FontSize', 12);
end
title('Drift-corrected model coefficients')
if saveFig
    pEnterStr = num2str(rm.modelParams.pEnter);
    saveFileName = ['driftCorrected_model_coeff_summary_stepThresh_', pEnterStr(3:end), ...
            '_', fileNameStr, '_interactions'];
    save_figure(f, figDir, saveFileName);
end
catch ME; rethrow(ME); end

%% Plot predicted vs. measured fluorescence
try
predictorVars = {'odorHistory', 'odorResp', 'fwSpeed'}; % odorHistory, odorResp, fwSpeed, yawSpeed
legendLocs = {'sw', 'sw', 'w', 'w'};
xLimits = {[], [], [300 600], [1800 3200], [700 1900], [300 1000], [300 600], [100 1100], [1000 1600]}';

legendStr = {'predicted fl', 'measured fl'};
for iExp = 1:size(rm.modelData, 1)
    
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
    xL = xLimits{iExp};
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
        suptitle([currModelData.expID{:}, '  —  drift-corrected model (adj. R^2 = ', ...
            num2str(currModelData.driftCorrectedMdlAdjR2, 2), ')'])
end

catch ME; rethrow(ME); end

% ==================================================================================================
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

%% Collect info from different model sets
try
fileIDStrings = {'01_with', '02_with', '03_with', '015_with', '025_with', '01_no', '02_no'};

dcR2 = [];
nohR2 = [];
dcNumCoeffs = [];
nohNumCoeffs = [];
bestWinSizes = [];
coeffNames = {};
coeffs = {};
for iFile = 1:numel(fileIDStrings)
    currFileName = fullfile(parentDir, ['rm_stepThresh_', fileIDStrings{iFile}, ...
            '_interactions.mat'])
      
    % Load dataset
    load(currFileName, 'rm');
    
    dcR2(iFile, :) = [rm.modelData.driftCorrectedMdlAdjR2'];
    nohR2(iFile, :) = [rm.modelData.noOdorHistMdlAdjR2'];
    dcNumCoeffs(iFile, :) = cellfun(@(x) x.NumCoefficients, rm.modelData.driftCorrectedMdls)';
    nohNumCoeffs(iFile, :) = cellfun(@(x) x.NumCoefficients, rm.modelData.noOdorHistMdls)';
    bestWinSizes(iFile, :) = rm.modelData.bestWinSizes';
    coeffNames{iFile} = cellfun(@(x) x.CoefficientNames, rm.modelData.driftCorrectedMdls, ...
            'UniformOutput', false);
    
end
catch ME; rethrow(ME); end