parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
dataDir = fullfile(parentDir, 'all_experiments');
lastExpID = []; lastOdorRespFilterDur = [];

%% Load and pre-process data for a new experiment or odorRespFilter duration

%Set data parameters
expID = '20190411-2';
roiName = 'TypeD';
trialNums = [];
maxSpeed = 100;
smWinVols = 3;
smWinFrames = 3;
smReps = 10;
ftLagVols = 3;
odorRespFilterDur = [7 5];

% Load and pre-process data
try
disp('Loading and preprocessing data...')
% -------------------- Load new data if necessary -------------------------%
if isempty(lastExpID) || ~strcmp(lastExpID, expID)
    % Load exp and trial metadata
    [expMd, trialMd] = load_metadata({expID}, dataDir);
    
    % Load imaging data
    roiData = load_roi_data({expID}, dataDir);
    roiNames = regexp(roiData.roiName, [roiName, '-.'], 'match');
    roiNames(cellfun(@isempty, roiNames)) = [];
    roiData = roiData(strcmp(roiData.roiName, roiNames{1}), :);
    
    % Load FicTrac data
    FRAME_RATE = 25;
    ft = load_ft_data({expID}, dataDir);
    
    % Load odor event data
    eventData = load_event_data({expID}, dataDir);
    odorEventData = eventData.odor;
    
    lastExpID = expID;
end
if ~exist('alignObj', 'var')
    % Load new alignEvent object
    disp('Loading new AlignEvent object...');
    load(fullfile(parentDir, 'Saved_AlignEvent_objects', '20200609_AlignEventObj_odor.mat'));
    disp('Loading complete')
end
if isempty(lastOdorRespFilterDur) || any(lastOdorRespFilterDur ~= odorRespFilterDur)
    % Load new event data table
    disp('Loading new event DataTable...')
    load(fullfile(parentDir, 'Saved_DataTable_objects', ['odor_pre_', num2str(odorRespFilterDur(1)), ...
            '_post_', num2str(odorRespFilterDur(2)), '.mat']), 'dt')    
    lastOdorRespFilterDur = odorRespFilterDur;
    disp('Loading complete');
end

% Process each trial (assuming no time between trials)
if isempty(trialNums)
   trialNums = 1:size(trialMd, 1); 
end
allFwSpeed = []; allYawSpeed = []; allOdorVols = []; allFl = []; allTrialStartDistances = [];
for iTrial = 1:numel(trialNums)
    currTrialNum = trialNums(iTrial);
    md = innerjoin(expMd, trialMd(currTrialNum, :));
    volTimes = calc_volTimes(md.nVolumes, md.volumeRate, md.trialDuration, md.originalTrialCount);
    
    % Smooth imaging rawFl data for current trial
    currRoiData = roiData(currTrialNum, :);
    currFl = smoothdata(currRoiData.rawFl{:}, 'gaussian', smWinVols);
    if ~isnan(md.pmtShutoffVols{:})
        currFl(md.pmtShutoffVols{:}) = nan;
%         currFl(1:max(md.pmtShutoffVols{:})) = nan; % Use to manually skip first part of experiment
    end
    
    % Smooth FicTrac data, then downsample to match volume rate
    currFt = ft(currTrialNum, :);
    currFwSpeed = repeat_smooth(currFt.fwSpeed{:} .* 4.5 .* FRAME_RATE, smReps, 'smWin', ...
            smWinFrames);
    currFwSpeed(currFwSpeed > maxSpeed) = maxSpeed;
    currYawSpeed = abs(repeat_smooth(rad2deg(currFt.yawSpeed{:}) .* FRAME_RATE, smReps, 'smWin', ...
            smWinFrames));
    currFwSpeedVols = []; currYawSpeedVols = [];
    for iVol = 1:numel(volTimes)
        dsFrame = argmin(abs(currFt.frameTimes{:} - volTimes(iVol)));
        currFwSpeedVols(iVol) = currFwSpeed(dsFrame);
        currYawSpeedVols(iVol) = currYawSpeed(dsFrame);
    end
    
    % Apply a lag to the FicTrac data
    currFwSpeedVols = currFwSpeedVols([ones(1, ftLagVols), 1:end-ftLagVols]);
    currYawSpeedVols = currYawSpeedVols([ones(1, ftLagVols), 1:end-ftLagVols]);
    
    % Get odor command vector for current trial
    currTrialOdorVols = odorEventData.create_logical_array(md.nVolumes, volTimes, table({expID}, ...
            currTrialNum, 'variableNames', {'expID', 'trialNum'}));
        
    % Append variables to whole-experiment vectors
    allFl = [allFl; currFl];
    allFwSpeed = [allFwSpeed; currFwSpeedVols'];
    allYawSpeed = [allYawSpeed; currYawSpeedVols'];
    allOdorVols = [allOdorVols; currTrialOdorVols];
    allTrialStartDistances = [allTrialStartDistances; (1:md.nVolumes)'];
    
end
catch ME; rethrow(ME); end

% Calculate the odor response filter from trial-averaged data
try
disp('Generating odor response vector...');
filterDefs = alignObj.create_filterDefs;
filterDefs.expID = expID;
filterDefs.roiName = roiName;
filterDefs.odor = -1;
filterDefs.locomotion = 0;
filterDefs.flailing = 0;
filterDefs.grooming = 0;
dt = dt.initialize_filters(filterDefs);
if size(dt.subset, 1) == 0
    % Relax the locomotion exclusion requirement if there isn't any data for this experiment
    filterDefs.locomotion = [];
    dt = dt.initialize_filters(filterDefs);
end
meanTrace = smoothdata(mean(cell2mat(dt.subset.eventFl'), 2, 'omitnan'), 'gaussian', smWinVols);
odorRespFilter = meanTrace(dt.subset.volTimes{1} > 0);
odorRespFilter = odorRespFilter - odorRespFilter(1);

% Plot filter shape
f = figure(22);clf; hold on;
f.Color = [1 1 1];
plot(dt.subset.volTimes{1}(dt.subset.volTimes{1} > 0)', odorRespFilter, 'linewidth', 3, 'color', ...
        'k');
plot_stim_shading([0, dt.subset.offsetTime(1) - dt.subset.onsetTime(1)])
title([regexprep(expID, '_', '\\_'), ' odor response filter'])
xlabel('Time (sec)');
set(gca, 'FontSize', 11)

% Extract onset volumes from odor command vector
odorOnsetVols = (regexp(regexprep(num2str(allOdorVols'), ' ', ''), '01')) + 1;
odorOnsetVector = zeros(size(allOdorVols));
odorOnsetVector(odorOnsetVols) = 1;

% Convolve with odor response filter to get predicted fl responses
odorRespVector = conv(odorOnsetVector, odorRespFilter);
odorRespVector = odorRespVector(1:numel(odorOnsetVector));
disp('Preprocessing complete');

catch ME; rethrow(ME); end

%% Generate and fit model

odorIntegrationWin = [26 28 30 32];
% odorIntegrationWin = [10 20 30 60 90 120];

% Model parameters
trainTestSplit = 0.8;
kFold = 100;
nSteps = [];
criterion = 'adjrsquared'; % 'sse, 'aic', 'bic', '', or 'adjrsquared'
pEnter = [0.01];
pRemove = [0];
verbose = 0;

% allFl(1:1500) = nan; % Manually throw out some of the data if the PMT gain was changed mid-exp

% Create table of predictor and response variables
tbl = table(abs(allFwSpeed), 'VariableNames', {'fwSpeed'});
tbl.yawSpeed = allYawSpeed;
tbl.odorResp = odorRespVector;
% tbl.volsFromTrialStart = allTrialStartDistances;
tbl.volsFromExpStart = (1:size(tbl, 1))';

% Prepare and split data
try   
% Calculate integrated odor values
odorIntegrationWinVols = round(odorIntegrationWin * md.volumeRate);
integratedOdor = [];
for iWin = 1:numel(odorIntegrationWinVols)
    currIntegrationWin = odorIntegrationWinVols(iWin);
    for iVol = 1:numel(allOdorVols)
        if iVol <= currIntegrationWin
            integratedOdor(iWin, iVol) = sum(allOdorVols(1:iVol));
        else
            integratedOdor(iWin, iVol) = sum(allOdorVols(iVol - currIntegrationWin:iVol));
        end
    end
end
for iWin = 1:size(integratedOdor, 1) 
    varName = ['odorHistory_', num2str(odorIntegrationWin(iWin))];
    tbl.(varName) = integratedOdor(iWin, :)';
end    

% Add fluorescence data in last column of input table
tbl.fl = allFl;
    
% Scale all variables so that the max value is 1;    
for iCol = 1:(size(tbl, 2))
    currData = tbl{:, iCol};
    tbl{:, iCol} = tbl{:, iCol} - mean(tbl{:, iCol}, 'omitnan'); % Center mean of all variables at zero
    tbl{:, iCol} = tbl{:, iCol} ./ max(tbl{:, iCol}, [], 'omitnan');
end

% % Calculate trial/chunk splitting indices
% splitInd = 1:size(tbl, 1);
% if md.originalTrialCount > 1
%     rsSplitInd = reshape(splitInd, size(tbl, 1) / sum(trialMd(trialNums, :).originalTrialCount), []);
% else
%     rsSplitInd = reshape(splitInd, numel(trialNums), []);
% end
% Split data into chunks/trials and add every other trial to test vs. predict datasets
% fitInds = as_vector(rsSplitInd(:, 1:2:end));
% predictInds = as_vector(rsSplitInd(:, 2:2:end));
% tblFit = tbl(fitInds, :);
% tblPredict = tbl(predictInds, :);

% Divide data into train/test groups
tblShuffle = tbl(randperm(size(tbl, 1)), :);
tblShuffle = tblShuffle(~isnan(tblShuffle.fl), :); % Drop rows with invalid fluorescence values
splitInd = floor(size(tblShuffle, 1) * trainTestSplit);
tblFit = tblShuffle(1:splitInd, :);
tblTest = tblShuffle((splitInd + 1):end, :);

% Divide training data into cross-validation folds
chunkSize = round(size(tblFit, 1) / kFold);
startInds = 1:chunkSize:size(tblFit, 1);
endInds = [startInds(2:end) - 1, size(tblFit, 1)];

catch ME; rethrow(ME); end

% --------------------- Model fitting ------------------------

% Select odor history window size using cross-validation
try
allPredAdjR2 = []; allMdlAdjR2 = []; allMdlAIC = []; allMdlBIC = []; allMdls = [];
for iWin = 1:numel(odorIntegrationWin)
    
    % Choose which odor integration window to use for current run
    odorHistoryCols = (size(tblFit, 2) - numel(odorIntegrationWin)):(size(tblFit, 2) - 1);
    currTbl = tblFit(:, [1:(odorHistoryCols(1) - 1), odorHistoryCols(iWin), size(tblFit, 2)]);
    
    % Fit model with cross-validation
    for iFold = 1:kFold
        if ~mod(iFold, 10)
            disp(['Running fold ', num2str(iFold), ' for window size = ', ...
                num2str(odorIntegrationWin(iWin))]);
        end 
        
        % Split data for current fold
        tblTrain = currTbl;
        tblTrain(startInds(iFold):endInds(iFold), :) = [];
        tblFoldTest = currTbl(startInds(iFold):endInds(iFold), :);
        
        % ------------- Generate model ---------------
        kvArgs = {'nSteps', nSteps, 'criterion', criterion, 'pEnter', pEnter, 'pRemove', pRemove, ...
                'verbose', verbose};
        emptyArgs = cellfun(@isempty, kvArgs);
        kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
        mdl = stepwiselm(tblTrain, kvArgs{:});
        % mdl = fitlm(tblFit);
        varNames = tblTrain.Properties.VariableNames();
        allMdls{iWin, iFold} = compact(mdl);
        
        % Evaluate model fit
        allMdlAdjR2(iWin, iFold) = mdl.Rsquared.Adjusted;
        allMdlAIC(iWin, iFold) = allMdls{iWin, iFold}.ModelCriterion.AIC;
        allMdlBIC(iWin, iFold) = allMdls{iWin, iFold}.ModelCriterion.BIC;
        
        % Evaluate fit to test data
        [yPred, ~] = predict(mdl, tblFoldTest(:, 1:end-1));
        [~, predAdjR2] = r_squared(tblFoldTest.fl, yPred, mdl.NumCoefficients);
        allPredAdjR2(iWin, iFold) = predAdjR2;
        
    end
    disp(' ');
    disp(' ');
    disp('--------------------------------------');
    disp(varNames{end - 1})
    disp('--------------------------------------');
    disp('Coefficient names:')
    namesMatch = 1;
    for iMdl = 1:kFold
        if iMdl == 1
            coeffNames = allMdls{iWin, iMdl}.CoefficientNames;
        else
            if numel(coeffNames) ~= numel(allMdls{iWin, iMdl}.CoefficientNames) || ...
                        ~all(strcmp(coeffNames, allMdls{iWin, iMdl}.CoefficientNames))
                namesMatch = 0;
            end
        end
    end
    if namesMatch
        disp(coeffNames);
    else
        disp('Warning: variable name mismatch across folds!');
    end
    disp(['Mean fit data adj R2: ', num2str(mean(allMdlAdjR2(iWin, :)), 2), ...
            ' (SD: ', num2str(std(allMdlAdjR2(iWin, :)), 2), ')']);
    disp(['Mean test data adj R2: ', num2str(mean(allPredAdjR2(iWin, :)), 2), ...
            ' (SD: ', num2str(std(allPredAdjR2(iWin, :)), 2), ')']);
    disp(['Mean model AIC: ', num2str(mean(round(allMdlAIC(iWin, :), 2))), ...
            ' (SD: ', num2str(std(allMdlAIC(iWin, :)), 2), ')']);
    disp(['Mean model BIC: ', num2str(mean(round(allMdlBIC(iWin, :), 2))), ...
            ' (SD: ', num2str(std(allMdlBIC(iWin, :)), 2), ')']);
    
end%iWin
disp('Testing complete');

% Summarize results
disp(' ')
disp('Fit R2:')
disp(num2str(mean(allMdlAdjR2, 2)', 2)); 
disp('Test R2:')
disp(num2str(mean(allPredAdjR2, 2)', 2)); 
disp('AIC:')
disp(num2str(round(mean(allMdlAIC, 2))'))
disp('Best odorIntegrationWin size by model adj. R2:')
disp(odorIntegrationWin(argmax(mean(allMdlAdjR2, 2), 1)));
disp('Best odorIntegrationWin size by test adj. R2:')
disp(odorIntegrationWin(argmax(mean(allPredAdjR2, 2), 1)));
disp('Best odorIntegrationWin size by AIC:')
disp(odorIntegrationWin(argmin(mean(allMdlAIC, 2), 1)));

catch ME; rethrow(ME); end

% Use the best window size to train and evaluate the final model
try

bestWinSize = argmax(mean(allPredAdjR2, 2), 1);
finalTrainTbl = tblFit(:, [1:(odorHistoryCols(1) - 1), odorHistoryCols(bestWinSize), ...
        size(tblFit, 2)]);
finalTestTbl = tblTest(:, [1:(odorHistoryCols(1) - 1), odorHistoryCols(bestWinSize), ...
        size(tblTest, 2)]);
finalMdl = stepwiselm(finalTrainTbl, kvArgs{:});
[yPred, ~] = predict(finalMdl, finalTestTbl(:, 1:end-1));
[~, predR2] = r_squared(finalTestTbl.fl, yPred, finalMdl.NumCoefficients);
predResiduals = finalTestTbl.fl - yPred;
predRMSE = mean(predResiduals.^2)^0.5;
disp(['Final model fit adj. R2: ', num2str(finalMdl.Rsquared.Adjusted, 2)])
disp(['Final model pred adj. R2: ', num2str(predR2, 2)])
disp(['Final model AIC: ', num2str(round(finalMdl.ModelCriterion.AIC))])
catch ME; rethrow(ME); end

%% --------------------- Plotting results ---------------------

% Plot coefficients for the final model
try
f = figure(12);clf; hold on;
f.Color = [1 1 1];
nCoeffs = finalMdl.NumCoefficients;
plot(2:nCoeffs, finalMdl.Coefficients.Estimate(2:end), 'o', 'color', 'k', 'linewidth', 3, 'markerSize', 8);
plot([0 nCoeffs + 1], [0 0], '--', 'color', 'b');
for iCoeff = 2:nCoeffs
    plot([iCoeff, iCoeff], [0, finalMdl.Coefficients.Estimate(iCoeff)], '-.', 'color', 'k')
end
ax = gca();
ax.XLim = [1 nCoeffs];
ax.FontSize = 14;
ax.XTick = 2:numel(finalMdl.Coefficients.Estimate);
ax.XTickLabel = regexprep(finalMdl.CoefficientNames(2:end), '_', '\\_');
ax.XTickLabelRotation = 45;
title([regexprep(expID, '_', '\\_'), ' model coefficients'])
catch ME; rethrow(ME); end

% Plot summary of CV training results
try
f = figure(10);clf;
f.Color = [1 1 1];
subaxis(1, 3, 1, 'mb', 0.15, 'ml', 0.06, 'mr', 0.03, 'sh', 0.08); hold on
set(gca, 'FontSize', 12)
plot(odorIntegrationWin, allMdlAdjR2, '-o');
plot(odorIntegrationWin, mean(allMdlAdjR2, 2), '-o', 'color', 'k', 'linewidth', 3)
title('Model fit adjusted R^2')
xlabel('odorInt window size')
ylabel('Adjusted R^2');
subaxis(1, 3, 2); hold on;
set(gca, 'FontSize', 12)
plot(odorIntegrationWin, allPredAdjR2, '-o');
plot(odorIntegrationWin, mean(allPredAdjR2, 2), '-o', 'color', 'k', 'linewidth', 3)
title('Model test data prediction adj. R^2')
xlabel('odorInt window size')
ylabel('Adjusted R^2');
subaxis(1, 3, 3); hold on;
set(gca, 'FontSize', 12)
plot(odorIntegrationWin, allMdlAIC, '-o');
plot(odorIntegrationWin, mean(allMdlAIC, 2), '-o', 'color', 'k', 'linewidth', 3)
title('Model AIC score')
xlabel('odorInt window size')
ylabel('AIC');
% subaxis(2, 3, 4)
% plot(allMdlAdjR2', '-o');
% xlabel('CV fold number')
% ylabel('Adjusted R2');
% subaxis(2, 3, 5)
% plot(allPredAdjR2', '-o');
% xlabel('CV fold number')
% ylabel('Adjusted R2');
% subaxis(2, 3, 6)
% plot(allMdlAIC', '-o');
% xlabel('CV fold number')
% ylabel('AIC');
catch ME; rethrow(ME); end

% Use final model to plot predicted fluorescence for all data 
predictorVars = {'odorHistory'}; % odorHistory, odorResp, fwSpeed, yawSpeed
try
f = figure(101); clf; hold on; 
f.Color = [1 1 1];

legendStr = {'predicted fl', 'measured fl'};
tblPred = tbl(:, [1:(odorHistoryCols(1) - 1), odorHistoryCols(bestWinSize), size(tbl, 2)]);
tblPred = tblPred(~isnan(tblPred.fl), :);
subaxis(1, 1, 1, 'mb', 0.14, 'mr', 0.05, 'ml', 0.05); hold on;
title([expID, '  —  Predicted vs. measured fluorescence (R^2 = ', num2str(predR2, 2), ')'])
plot(volTimes(1:size(tblPred, 1)), predict(finalMdl, tblPred(:, 1:end-1)), ...
        'linewidth', 1, 'color', rgb('blue'));
plot(volTimes(1:size(tblPred, 1)), tblPred.fl, 'linewidth', 1, 'color', rgb('orangered'));
yyaxis right; hold on
if any(strcmpi('odorResp', predictorVars))
   plot(volTimes(1:size(tblPred, 1)), tblPred.odorResp, 'color', rgb('darkmagenta'));
   legendStr{end + 1} = 'Odor response';
end
if any(strcmpi('fwSpeed', predictorVars))
   plot(volTimes(1:size(tblPred, 1)), tblPred.fwSpeed, '-', 'color', rgb('green'));
   legendStr{end + 1} = 'abs(fwSpeed)';
end
if any(strcmpi('yawSpeed', predictorVars))
   plot(volTimes(1:size(tblPred, 1)), tblPred.yawSpeed, 'color', rgb('red'));
   legendStr{end + 1} = 'yawSpeed';
end
if any(strcmpi('odorHistory', predictorVars))
    varName = ['odorHistory_', num2str(odorIntegrationWin(bestWinSize))];
    plot(volTimes(1:size(tblPred, 1)), tblPred.(varName), '-', 'color', rgb('gold'));
    legendStr{end + 1} = regexprep(varName, '_', '\\_');
end
legend(legendStr, 'fontsize', 14, 'Location', 'best')
ax = gca();
ax.FontSize = 14;
xlabel('Time (sec)')
catch ME; rethrow(ME); end

% Show scatter plots of predictor variables vs. fluorescence
try
f = figure(102); clf;
f.Color = [1 1 1];
varNames = finalMdl.PredictorNames;
nPlots = numSubplots(numel(varNames));
tblPred = tbl(:, [1:(odorHistoryCols(1) - 1), odorHistoryCols(bestWinSize), size(tbl, 2)]);
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
    plt = plot(tbl.(varNames{iPlot}), fit, '-', 'color', 'r', 'linewidth', 2.5);
    R = corrcoef(tblPred.(varNames{iPlot}), tblPred.fl);
    Rsquared = R(1,2)^2;
    legend(plt, [' y = ', num2str(round(p(1), 2)), 'x,  R^2 = ', num2str(Rsquared, 2)], ...
            'fontSize', 12)
    xlabel(regexprep(varNames{iPlot}, '_', '\\_'));
    ylabel('Fluorescence');
    ax.FontSize = 14;
end
catch ME; rethrow(ME); end
