parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
dataDir = fullfile(parentDir, 'all_experiments');
 lastOdorRespFilterDur = [];
 
%%

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
                     
          
expInfoTbl = table(expIDList', skipTrials', skipVols', 'VariableNames', {'expID', 'skipTrials', ...
        'skipVols'});

% Create analysis object
rm = RegressionModelAnalysis(expInfoTbl, p);

% Initialize models
rm = rm.initialize_models(mp);

% Find optimal odor integration windows for each experiment
rm = rm.optimize_integration_windows();

%% Use the calculated best window sizes to fit new models for each of the experiments 
fullMdls = {};
fullMdlPredFl = {};
fullMdlPredAdjR2 = [];
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
    bestHistVarInd = ~cellfun(@isempty, regexp(varNames, ['_', num2str(bestWinSize)]));
    tblFit = currModelData.fullDataTbl{:}(currModelData.fitRowInds{:}, ...
            logical(~odorHistVarInds + bestHistVarInd));
    tblTest = currModelData.fullDataTbl{:}(currModelData.testRowInds{:}, logical(~odorHistVarInds ...
            + bestHistVarInd));

    % Use all training data to create and evaluate a final stepwise model
    kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
            mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
    emptyArgs = cellfun(@isempty, kvArgs);
    kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
    fullMdls{iExp} = stepwiselm(tblFit, kvArgs{:});
    [fullMdlPredFl{iExp}, ~] = predict(fullMdls{iExp}, tblTest(~isnan(tblTest.fl), 1:end-1));
    [~, fullMdlPredAdjR2(iExp)] = r_squared(tblTest.fl(~isnan(tblTest.fl)), fullMdlPredFl{iExp}, ...
            fullMdls{iExp}.NumCoefficients);
end
disp('Final models created');
rm.modelData.fullMdls = fullMdls';
rm.modelData.finalMdlAdjR2 = fullMdlPredAdjR2';
disp(rm.modelData)

%% Plot model coefficients

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

%% Subtract volsFromOdorStart from fl data and re-fit models

finalMdls = {};
predAdjR2 = [];
predFl = {};
disp('Training final models...')
allVolAdjTblTrain = {};
allVolAdjTblTest = {};
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
    bestHistVarInd = ~cellfun(@isempty, regexp(varNames, ['_', num2str(bestWinSize)]));
    tblFit = currModelData.fullDataTbl{:}(currModelData.fitRowInds{:}, ...
        logical(~odorHistVarInds + bestHistVarInd));
    tblTest = currModelData.fullDataTbl{:}(currModelData.testRowInds{:}, logical(~odorHistVarInds ...
        + bestHistVarInd));
    
    currMdl = rm.modelData.finalMdls{iExp};
    if ismember('volsFromExpStart', currMdl.Coefficients.Properties.RowNames)
        tblFitVolAdjust = tblFit;
        tblFitVolAdjust.fl = tblFitVolAdjust.fl - (currMdl.Coefficients{'volsFromExpStart', ...
            'Estimate'} .* tblFitVolAdjust.volsFromExpStart);
        tblTestVolAdjust = tblTest;
        tblTestVolAdjust.fl = tblTestVolAdjust.fl - (currMdl.Coefficients{'volsFromExpStart', ...
            'Estimate'} .* tblTestVolAdjust.volsFromExpStart);
        
        
        tblFit = tblFitVolAdjust(:, ~strcmp(tblFitVolAdjust.Properties.VariableNames, 'volsFromExpStart'));
        tblTest = tblTestVolAdjust(:, ~strcmp(tblTestVolAdjust.Properties.VariableNames, 'volsFromExpStart'));
    end    
    
    % Use all training data to create and evaluate a final stepwise model
    kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
        mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
    emptyArgs = cellfun(@isempty, kvArgs);
    kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
    finalMdls{iExp} = stepwiselm(tblFit, kvArgs{:});
    [yPred, ~] = predict(finalMdls{iExp}, tblTest(~isnan(tblTest.fl), 1:end-1));
    [~, predAdjR2(iExp)] = r_squared(tblTest.fl(~isnan(tblTest.fl)), yPred, ...
            finalMdls{iExp}.NumCoefficients);
    predFl = yPred;
end
disp('Final models created');
% rm.modelData.finalMdls = finalMdls';
% rm.modelData.finalMdlAdjR2 = predAdjR2';
% disp(rm.modelData)

% Plot coefficients
f = figure(2); clf; hold on;
f.Color = [1 1 1];
subplotDims = [2 5];
for iExp = 1:numel(finalMdls)
    
    currMdl = finalMdls{iExp};
    
    subaxis(subplotDims(1), subplotDims(2), iExp, 'mt', 0.05, 'mb', 0.12, 'sv', 0.16, ...
            'ml', 0.05, 'mr', 0.05);
    hold on
    ax = gca();
    RegressionModelAnalysis.plot_coeffs(currMdl, ax);
    title([regexprep(rm.modelData.expID{iExp}, '_', '\\_'), ' coeffs (adj. R^2 = ', ...
            num2str(predAdjR2(iExp), 2), ')'])
end

%% Now fit the "drift-corrected" models without the odor history variables


noHistMdls = {};
noHistPredAdjR2 = [];
noHistPredFl = {};
disp('Training no-odor history models...')
allVolAdjTblTrain = {};
allVolAdjTblTest = {};
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
    tblFit = currModelData.fullDataTbl{:}(currModelData.fitRowInds{:}, ~odorHistVarInds);
    tblTest = currModelData.fullDataTbl{:}(currModelData.testRowInds{:}, ~odorHistVarInds);
    
    % Subtract volsFromExpStart offset from fluorescence data
    currMdl = rm.modelData.finalMdls{iExp};
    if ismember('volsFromExpStart', currMdl.Coefficients.Properties.RowNames)
        tblFitVolAdjust = tblFit;
        tblFitVolAdjust.fl = tblFitVolAdjust.fl - (currMdl.Coefficients{'volsFromExpStart', ...
            'Estimate'} .* tblFitVolAdjust.volsFromExpStart);
        tblTestVolAdjust = tblTest;
        tblTestVolAdjust.fl = tblTestVolAdjust.fl - (currMdl.Coefficients{'volsFromExpStart', ...
            'Estimate'} .* tblTestVolAdjust.volsFromExpStart);
        
        
        tblFit = tblFitVolAdjust(:, ~strcmp(tblFitVolAdjust.Properties.VariableNames, 'volsFromExpStart'));
        tblTest = tblTestVolAdjust(:, ~strcmp(tblTestVolAdjust.Properties.VariableNames, 'volsFromExpStart'));
    end
    
    % Use all training data to create and evaluate a final stepwise model
    kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
        mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
    emptyArgs = cellfun(@isempty, kvArgs);
    kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
    noHistMdls{iExp} = stepwiselm(tblFit, kvArgs{:});
    [yPred, ~] = predict(noHistMdls{iExp}, tblTest(~isnan(tblTest.fl), 1:end-1));
    [~, noHistPredAdjR2(iExp)] = r_squared(tblTest.fl(~isnan(tblTest.fl)), yPred, ...
            noHistMdls{iExp}.NumCoefficients);
    noHistPredFl = yPred;
end
disp('Final models created');
% rm.modelData.finalMdls = finalMdls';
% rm.modelData.finalMdlAdjR2 = predAdjR2';
% disp(rm.modelData)

% Plot coefficients
f = figure(3); clf; hold on;
f.Color = [1 1 1];
subplotDims = [2 5];
for iExp = 1:numel(noHistMdls)
    
    currMdl = noHistMdls{iExp};
    
    subaxis(subplotDims(1), subplotDims(2), iExp, 'mt', 0.05, 'mb', 0.12, 'sv', 0.16, ...
            'ml', 0.05, 'mr', 0.05);
    hold on
    ax = gca();
    RegressionModelAnalysis.plot_coeffs(currMdl, ax);
    title([regexprep(rm.modelData.expID{iExp}, '_', '\\_'), ' coeffs (adj. R^2 = ', ...
            num2str(noHistPredAdjR2(iExp), 2), ')'])
    
end



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
        [yPred, ~] = predict(mdl, currFoldTestTbl(:, 1:end-1));
        [~, R2] = r_squared(currFoldTestTbl.fl, yPred, mdl.NumCoefficients);
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
[yPred, ~] = predict(finalMdl, tblTest(~isnan(tblTest.fl), 1:end-1));
[~, predAdjR2] = r_squared(tblTest.fl(~isnan(tblTest.fl)), yPred, ...
        finalMdl.NumCoefficients);

disp('Final model created');


catch ME; rethrow(ME); end


