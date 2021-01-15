parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
% dataDir = fullfile(parentDir, 'all_experiments');
dataDir = fullfile(parentDir, 'new_PPL201_experiments');
figDir = ...
        'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Figs\linear_regression_analysis';

%% Load one or more existing model objects
threshStr = {'03'};%{'01', '015', '02'};
allRm = {};
for i=1:numel(threshStr)
    load(fullfile(parentDir, 'Saved_LinearModel_objects', ...
            ['rm_stepThresh_', threshStr{i}, '_with_interactions_30_sec_min_winSize.mat']), 'rm');
    allRm{i} = rm;
end

%% Set up a new model 

% Set source data parameters;
p = [];
p.roiName = 'TypeB';
p.maxSpeed = 100;
p.smWinVols = 5;
p.smWinFrames = 5;
p.smReps = 15;
p.ftLagVols = 3;
p.speedType = 'moveSpeed';
p.odorRespFilterDur = [6 7];

p.parentDir = parentDir;
p.dataDir = dataDir;

p.eventDataParentDir = p.dataDir;
p.alignEventDateStr = '20210104';

% p.eventDataParentDir = p.parentDir;
% p.alignEventDateStr = '20200609';

p.convertFtUnits = 0;
p.loadOneNoteData = 0;
p.alignObjFilterDefs = [];

% Set model parameters
mp = [];
mp.trainTestSplit = 0.8;
mp.kFold = 100;
mp.criterion = 'rsquared'; % 'sse, 'aic', 'bic', 'rsquared', or 'adjrsquared'
mp.upper = [];
mp.pEnter = [0];
mp.pRemove = [0];
mp.verbose = 0;
mp.useYaw = 1;
mp.useDriftCorrection = 1;
% mp.odorIntegrationWin = [30:20:200];
mp.odorIntegrationWin = [30:60:600];
mp.speedPadDist = 5;
mp.speedIntegrationWin = [];
mp.standardizeInputs = 1;
mp.normalizeInputs = 0;


% % PPL203
% expIDList = {'20190304-1', '20190315-1', '20190315-2', '20190315-3', '20190401-1', ... 
%                '20190401-2', '20190411-1', '20190411-2', '20190411-3'};       
% skipTrials = {[], [], [], [], [], [], [], [], []};
% skipVols = {[], [1:1500], [], [], [], [], [], [], [1:2200]};

            
% % PPM1201
% expIDList = {'20180329-2', '20180405-2', '20180414-1', '20180414-2', '20180416-1', '20180523-2', ...
%         '20181020-1', '20190226-3'};
% skipTrials = repmat({[]}, 1, 8);
% skipVols = skipTrials;

% PPL201
expIDList = {'20201222-1', '20201222-2', '20201228-1', '20201228-2', '20201228-3', '20210102-1', ...
        '20210102-2', '20210102-3', '20210102-4'};
skipTrials = {[6], [6], [6], [6], ...
              [4:6], [6], [6], [6 7], [6 7]};      
% skipTrials = {[1 3:6], [1 3:6], [1 3:6], [1 3:6], ...
%               [1 3:6], [1 3:6], [1 3:6], [1 3:6 7], [1 3:6 7]};
skipVols = repmat({[]}, 1, numel(skipTrials));




                     
try          
expInfoTbl = table(expIDList', skipTrials', skipVols', 'VariableNames', {'expID', 'skipTrials', ...
        'skipVols'});

% Create analysis objects
baseRm = RegressionModelAnalysis(expInfoTbl, p);
criterionList = {'adjrsquared'};
pEnterList = [0.01];
pRemoveList = [0];
allRm = [];
for i = 1:numel(criterionList)
    disp(['Fitting models using ', criterionList{i}, ' criterion']);
    mp.criterion = criterionList{i};
    mp.pEnter = pEnterList(i);
    mp.pRemove = pRemoveList(i);
    allRm{i} = baseRm.initialize_models(mp);
    allRm{i} = allRm{i}.optimize_odor_integration_windows();
end

% % Create analysis object
% rm = RegressionModelAnalysis_PPL201(expInfoTbl, p);
% 
% % Initialize models
% rm = rm.initialize_models(mp);
% 
% % Find optimal odor integration windows for each experiment
% rm = rm.optimize_odor_integration_windows();

% Find optimal fw speed integration window sizes
% rm = rm.optimize_mean_speed_windows();

% allRm_normalized = allRm;
% save(fullfile(parentDir, 'model_testing_normalized.mat'), 'baseRm', 'allRm_normalized');

catch ME; rethrow(ME); end
%% =================================================================================================
%  MODEL POSTPROCESSING
%  =================================================================================================

%% Use the calculated best odor window sizes to fit new models for each of the experiments 
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

                    %% Use best speed history window sizes instead of odor integration windows
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
    bestWinSize = 10;
    varNames = currModelData.fullDataTbl{:}.Properties.VariableNames;
    speedHistVarInds = ~cellfun(@isempty,regexp(varNames, 'SpeedHistory')); 
    bestHistVarInd = ~cellfun(@isempty, regexp(varNames, ['_', num2str(bestWinSize), '$']));
    tblFit = currModelData.fullDataTbl{:}(currModelData.fitRowInds{:}, ...
            logical(~speedHistVarInds + bestHistVarInd));
    tblTest = currModelData.fullDataTbl{:}(currModelData.testRowInds{:}, logical(~speedHistVarInds ...
            + bestHistVarInd));
    tblPred = currModelData.fullDataTbl{:}(:, logical(~speedHistVarInds + bestHistVarInd));

%     expStartColInd = strcmp(tblFit.Properties.VariableNames, 'volsFromExpStart');
%     tblFit = tblFit(:, ~expStartColInd);
%     tblTest = tblTest(:, ~expStartColInd);
%     tblPred = tblPred(:, ~expStartColInd);    
    
    % Use all training data to create and evaluate a final stepwise model
    kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
            mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
    emptyArgs = cellfun(@isempty, kvArgs);
    kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
    fullMdls{iExp} = stepwiselm(tblFit, kvArgs{:});
    
%     fullMdls{iExp} = fitlm(tblFit, 'fl ~ 1 + fwSpeed + fwSpeed:fwSpeedHistory_10');
    
    
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

noSpeedHistMdls = {};
noSpeedHistAdjR2 = {};
for iExp = 1:size(rm.modelData, 1)
        if numel(rm.modelParams) > 1
        mp = rm.modelParams(iExp);
    else
        mp = rm.modelParams;
    end
    disp(rm.sourceData.expID{iExp})
    
    % Select data with best odorIntegration win for current experiment
    currModelData = rm.modelData(iExp, :);
    bestWinSize = 120;
    varNames = currModelData.fullDataTbl{:}.Properties.VariableNames;
    speedHistVarInds = ~cellfun(@isempty,regexp(varNames, 'SpeedHistory')); 
    bestHistVarInd = ~cellfun(@isempty, regexp(varNames, ['_', num2str(bestWinSize), '$']));
    tblFit = currModelData.fullDataTbl{:}(currModelData.fitRowInds{:}, ...
            logical(~speedHistVarInds + bestHistVarInd));
    tblTest = currModelData.fullDataTbl{:}(currModelData.testRowInds{:}, logical(~speedHistVarInds ...
            + bestHistVarInd));
    tblPred = currModelData.fullDataTbl{:}(:, logical(~speedHistVarInds + bestHistVarInd));

    expStartColInd = strcmp(tblFit.Properties.VariableNames, 'volsFromExpStart');
%     tblFit = tblFit(:, ~expStartColInd);
%     tblTest = tblTest(:, ~expStartColInd);
%     tblPred = tblPred(:, ~expStartColInd);
         
    noSpeedHistMdls{iExp} = removeTerms(fullMdls{iExp}, 'fwSpeedHistory_10');
    
    [predFl, ~] = predict(noSpeedHistMdls{iExp}, tblTest(~logical(sum(isnan(table2array(tblTest)), 2)), ...
            1:end-1));
    [~, noSpeedHistMdlAdjR2(iExp)] = r_squared(tblTest.fl(~logical(sum(isnan(table2array(tblTest)), 2))), ...
            predFl, noSpeedHistMdls{iExp}.NumCoefficients);
%     [fullMdlPredFl{iExp}, ~] = predict(noSpeedHistMdls{iExp}, ...
%             tblPred(~logical(sum(isnan(table2array(tblPred)), 2)), 1:end-1));
end

rm.modelData.noSpeedHistMdls = noSpeedHistMdls';
rm.modelData.noSpeedHistMdlAdjR2 = noSpeedHistMdlAdjR2';



catch ME; rethrow(ME); end
                    %% Train initial models for speed history windows
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
    tblFit = currModelData.fullDataTbl{:}(currModelData.fitRowInds{:}, :);
    tblTest = currModelData.fullDataTbl{:}(currModelData.testRowInds{:}, :);
    tblPred = currModelData.fullDataTbl{:}(:, :); 
    
    % Use all training data to create and evaluate a stepwise model
    kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
            mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
    emptyArgs = cellfun(@isempty, kvArgs);
    kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
    fullMdls{iExp} = stepwiselm(tblFit, kvArgs{:});
    
%     fullMdls{iExp} = fitlm(tblFit, 'fl ~ 1 + fwSpeed + fwSpeed:fwSpeedHistory_10');

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
    driftCorrectedMdlMeasuredFl{iExp} = tblPred.fl(~logical(sum(isnan(table2array(tblPred)), 2)));
    driftCorrectedMdlTestFl{iExp} = tblTest.fl(~logical(sum(isnan(table2array(tblTest)), 2)));
    
    % Use all training data to create and evaluate a final stepwise model
    kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
        mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
    emptyArgs = cellfun(@isempty, kvArgs);
    kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
    driftCorrectedMdls{iExp} = stepwiselm(tblFit, kvArgs{:});
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
    [predFl, ~] = predict(noOdorHistReFitMdls{iExp}, ...
            tblTest(~logical(sum(isnan(table2array(tblTest)), 2)), 1:end-1));
    [~, noOdorHistReFitMdlAdjR2(iExp)] = r_squared(rm.modelData.driftCorrectedMdlTestFl{iExp}, predFl, ...
            noOdorHistReFitMdls{iExp}.NumCoefficients);
    [noOdorHistReFitMdlPredFl{iExp}, ~] = predict(noOdorHistReFitMdls{iExp}, ...
           tblPred(~logical(sum(isnan(table2array(tblPred)), 2)), 1:end-1));
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
    [predFl, ~] = predict(noOdorHistMdls{iExp}, ...
            tblTest(~logical(sum(isnan(table2array(tblTest)), 2)), ...
            ~ismember(tblTest.Properties.VariableNames, odorHistCoeffs)));
    [~, noOdorHistMdlAdjR2(iExp)] = r_squared( ...
            tblTest.fl(~logical(sum(isnan(table2array(tblTest)), 2))), predFl, ...
            noOdorHistMdls{iExp}.NumCoefficients);
    [noOdorHistMdlPredFl{iExp}, ~] = predict(noOdorHistMdls{iExp}, ...
           tblPred(~logical(sum(isnan(table2array(tblPred)), 2)), 1:end-1));
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

%% Plot overlaid odor response filters and trial-averaged behavioral responses to odor
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
plotArr = [rm.modelData.driftCorrectedMdlAdjR2, ...
        rm.modelData.noOdorHistMdlAdjR2]';
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
ylim([0, 1])
ax = gca();
ax.XTick = 1:size(plotArr, 1);
ax.XTickLabel = {'Drift-corrected', 'No odor history', 'No odor history (re-fit)'};
% ax.XTickLabel = {'With interactions', 'No interactions', 'No odor history'};
ylabel('Model adjusted R^2');
ax.FontSize = 14;
% titleStr = '';
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
saveFig = 0;
try
allDcCoeffNames = {};
for iExp = 1:size(rm.modelData, 1)
    allDcCoeffNames = [allDcCoeffNames, rm.modelData.driftCorrectedMdls{iExp}.CoefficientNames];    
end
uniqueCoeffs = unique(regexprep(allDcCoeffNames, 'odorHistory_.*', 'odorHistory'));

uniqueCoeffs = uniqueCoeffs(2:end) % Drop intercept term
uniqueCoeffs = uniqueCoeffs([1 4 3 5 2]);
% uniqueCoeffs = uniqueCoeffs([1 3 2]);
coeffArrDc = zeros(numel(uniqueCoeffs), size(rm.modelData, 1));
for iExp = 1:size(rm.modelData, 1) 
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
coeffTblDc = array2table(coeffArrDc', 'variableNames', uniqueCoeffs, 'rowNames', ...
        rm.modelData.expID);
f = figure(35); clf; 
f.Color = [1 1 1];
f.Position = [f.Position(1:2), 1035, 950];
if (f.Position(2) + f.Position(4)) > 1000
    f.Position(2) = 1000 - f.Position(4);
end
plotArr = table2array(coeffTblDc);
lim = max(abs([max(plotArr(:)), min(plotArr(:))]));
xSize = size(plotArr, 2);
plotArr(:, end + 1) = -lim;
plotArr(:, end + 1) = lim;
ax = subaxis(1, 1, 1, 'ml', 0.12, 'mt', 0.12, 'mb', 0.03, 'mr', 0.11);
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
    text(2.5, iRow, ['   Adj. R^2 = ', num2str(rm.modelData.driftCorrectedMdlAdjR2(iRow), 2)], ...
            'FontSize', 12);
end
title('Drift-corrected model coefficients')
f.UserData.sourceDataParams = rm.sourceDataParams;
f.UserData.modelParams = rm.modelParams;
if saveFig
    pEnterStr = num2str(rm.modelParams.pEnter);
    saveFileName = ['model_coeff_summary_driftCorrected_stepThresh_', pEnterStr(3:end), ...
            '_', fileNameStr, '_interactions'];
    save_figure(f, figDir, saveFileName);
end
catch ME; rethrow(ME); end

%% Plot predicted vs. measured fluorescence

saveFigs = 0;
fileNameStr = 'fullExp_';


predictorVars = {'odorHistory', 'odorResp', 'moveSpeed'}; % odorHistory, odorResp, fwSpeed, yawSpeed
% legendLocs = {'sw', 'sw', 'nw', 'nw'};
legendLocs = {'best', 'best', 'best', 'best'};
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
    if any(strcmpi('moveSpeed', predictorVars))
        plot(xx, fullDataTbl.moveSpeed(~logical(sum(isnan(table2array(fullDataTbl)), 2))), '-', ...
                'color', rgb('blue'), 'linewidth', 1);
        legendStr{end + 1} = 'abs(fwSpeed)';
    end
    if any(strcmpi('yawSpeed', predictorVars))
        plot(xx, fullDataTbl.yawSpeed(~logical(sum(isnan(table2array(fullDataTbl)), 2))), 'color', ...
                rgb('red'));
        legendStr{end + 1} = 'yawSpeed';
    end
    if any(strcmpi('odorHistory', predictorVars))
        varName = ['odorHistory_', num2str(bestWinSize)];
        plot(xx, fullDataTbl.(varName)(~logical(sum(isnan(table2array(fullDataTbl)), 2))), '-', ...
                'color', rgb('black'), 'linewidth', 3);
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
        
    if saveFigs
        figTitle = ['predictedVsMeasuredFl_', fileNameStr, currModelData.expID{:}];
        save_figure(f, figDir, figTitle);
    end
end

catch ME; rethrow(ME); end

% ==================================================================================================







