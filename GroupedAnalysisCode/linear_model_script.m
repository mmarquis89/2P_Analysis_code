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
odorRespFilterDur = [2 6];

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

% Generate odor response filter and vector
disp('Generating odor response vector...');
% ------------------- Convolution ------------------------%
% Create odor response filter
% gaussian = @(a, b, c, x) a .* exp(-( (x - b).^2 ./ (2 * c^2) ));
% a1 = 0.3;
% b1 = 1 .* md.volumeRate;
% c1 = 0.3 .* md.volumeRate;
% a2 = 1;
% b2 = 3 .* md.volumeRate;
% c2 = 1 .* md.volumeRate;
% a1 = 0.2;
% b1 = 0.6 .* md.volumeRate;
% c1 = 0.2 .* md.volumeRate;
% a2 = 1;
% b2 = 2.5 .* md.volumeRate;
% c2 = 0.8 .* md.volumeRate;
% g1 = gaussian(a1, b1, c1, 1:1:(round(6 * md.volumeRate)));
% g2 = gaussian(a2, b2, c2, 1:1:(round(6 * md.volumeRate)));

% figure(1);clf; hold on;
% plot(volTimes(1:(round(6 * md.volumeRate))), g1);
% plot(volTimes(1:(round(6 * md.volumeRate))), g2);
% plot(volTimes(1:(round(6 * md.volumeRate))), g1 + g2, 'linewidth', 3);
% plot_stim_shading([0,4])
% xlim([-2, 4])
% odorRespFilter = [0, g1 + g2];

% CALCULATE ODOR RESPONSE FILTER FROM TRIAL-AVERAGED DATA
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
figure(22);clf; hold on;
plot(dt.subset.volTimes{1}(dt.subset.volTimes{1} > 0)', odorRespFilter, 'linewidth', 3);
plot_stim_shading([0, dt.subset.offsetTime(1) - dt.subset.onsetTime(1)])

% Extract onset volumes from odor command vector
odorOnsetVols = (regexp(regexprep(num2str(allOdorVols'), ' ', ''), '01')) + 1;
odorOnsetVector = zeros(size(allOdorVols));
odorOnsetVector(odorOnsetVols) = 1;

% Convolve with odor response filter to get predicted fl responses
odorRespVector = conv(odorOnsetVector, odorRespFilter);
odorRespVector = odorRespVector(1:numel(odorOnsetVector));

catch ME; rethrow(ME); end

%% Generate and fit model


odorIntegrationWin = [30];

% Model parameters
trainSplit = 0.25;
nSteps = [];
criterion = 'adjrsquared'; % 'sse, 'aic', 'bic', 'rsquared', or 'adjrsquared'
pEnter = [0.01];
pRemove = [];
verbose = 2;




% Create table of predictor and response variables
tbl = table(allFwSpeed, 'VariableNames', {'fwSpeed'});
tbl.yawSpeed = allYawSpeed;
tbl.odorResp = odorRespVector;
% tbl.dummyVar = allFl(randperm(numel(allFl))); tbl.dummyVar(isnan(tbl.dummyVar)) = 0;
% tbl.volsFromTrialStart = allTrialStartDistances;
tbl.volsFromExpStart = (1:size(tbl, 1))';

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
% allFl(1:1800) = nan;

% Divide data into train/test groups
% try
    
% Calculate trial/chunk splitting indices
splitInd = 1:size(tbl, 1);
if md.originalTrialCount > 1
    rsSplitInd = reshape(splitInd, size(tbl, 1) / sum(trialMd(trialNums, :).originalTrialCount), []);
else
    rsSplitInd = reshape(splitInd, numel(trialNums), []);
end

% Scale all variables so that the max value is 1;    
for iCol = 1:(size(tbl, 2))
    currData = tbl{:, iCol};
    tbl{:, iCol} = tbl{:, iCol} - mean(tbl{:, iCol}, 'omitnan');% Center all variables at zero
    tbl{:, iCol} = tbl{:, iCol} ./ max(tbl{:, iCol}, [], 'omitnan');
end
     
% Split data into chunks/trials and add every other trial to test vs. predict datasets
% fitInds = as_vector(rsSplitInd(:, 1:2:end));
% predictInds = as_vector(rsSplitInd(:, 2:2:end));
% tblFit = tbl(fitInds, :);
% tblPredict = tbl(predictInds, :);

% Divide test and predict data by individual volumes instead of trials
tblShuffle = tbl(randperm(size(tbl, 1)), :);
tblShuffle = tblShuffle(~isnan(tblShuffle.fl), :);
splitPoint = floor(size(tblShuffle, 1) * trainSplit); % Drop rows with invalid fluorescence values
tblFit = tblShuffle(1:splitPoint, :);
tblPredict = tblShuffle((splitPoint + 1):end, :);

% catch ME; rethrow(ME); end

% ------------- Generate model ---------------
kvArgs = {'nSteps', nSteps, 'criterion', criterion, 'pEnter', pEnter, 'pRemove', pRemove, ...
        'verbose', verbose};
emptyArgs = cellfun(@isempty, kvArgs);
kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];

% Fit model
mdl = stepwiselm(tblFit, kvArgs{:});
% mdl = fitlm(tblFit);
varNames = tblFit.Properties.VariableNames();
disp(mdl)

%

% for i = 1:(numel(varNames) - 1)
%     figure(i); clf; hold on; plot(tbl.(varNames{i}), tbl.fl, '.');
%     p = polyfit(tbl.(varNames{i})(~isnan(tbl.fl)), tbl.fl(~isnan(tbl.fl)), 1);
%     f = polyval(p, tbl.(varNames{i}), 'omitnan'); plot(tbl.(varNames{i}), f, '-');
%     legend('Raw data', [num2str(round(p(1), 2)), 'x + ', num2str(round(p(2), 2))])
%     xlabel(varNames{i}); ylabel('fl'); title('Raw data scatter')
%     figure(i + numel(varNames) - 1); clf; plotAdded(mdl, varNames{i}); 
% %     figure(i + 2 * (numel(varNames) - 1)); clf; plotAdjustedResponse(mdl, varNames{i});
% end
figure(100);clf; plotEffects(mdl);

% Use model to predict fluorescence for the rest of the data
startInds = 1:size(tblFit, 1):size(tblPredict, 1);
allPredR2 = [];
allPredRMSE = [];
allPredAIC = [];
for iFold = 1:numel(startInds)
    currStart = startInds(iFold);
    if iFold < numel(startInds)
        currEnd = startInds(iFold + 1) - 1;
    else
        currEnd = size(tblPredict, 1);
    end
    currFoldPredict = tblPredict(currStart:currEnd, :);
    [yPred, ~] = predict(mdl, currFoldPredict(:, 1:end-1));
    [~, predR2] = r_squared(currFoldPredict.fl, yPred, mdl.NumCoefficients);
    predResiduals = currFoldPredict.fl - yPred;
    predRMSE = mean(predResiduals.^2)^0.5;
    allPredR2(iFold) = predR2;
    allPredRMSE(iFold) = predRMSE;
    allPredAIC(iFold) = aicbic(mdl.LogLikelihood, mdl.NumCoefficients);
end
disp(['Fit data adj R2: ', num2str(mdl.Rsquared.adjusted, 2)]);
disp(['Test data adj R2s: ', num2str(allPredR2, 2)]);
disp(['Mean = ', num2str(mean(allPredR2), 2)]);
% disp(['Fit data AIC: ', num2str(aicbic(mdl.LogLikelihood, mdl.NumCoefficients), 3)])
% disp(['Test data AICs: ', num2str(allPredAIC, 3)])
% disp(['Mean = ', num2str(mean(allPredAIC), 3)]);
% 
% fprintf(['\nFit data adj R2: ', num2str(mdl.Rsquared.adjusted, 3), ', RMSE: ', num2str(mdl.RMSE, 3), ...
%         '\nTest data adj R2: ', num2str(predR2, 3), ', RMSE: ', num2str(predRMSE, 3), ...
%         '\n']);
% [aic, bic] = aicbic(mdl.LogLikelihood, mdl.NumCoefficients, size(tblFit, 1));
% disp(['AIC: ', num2str(aic)])
% disp(['BIC: ', num2str(bic)])

figure(101);clf;hold on; 
plot(predict(mdl, tbl(~isnan(tbl.fl), 1:end-1)));
plot(tbl.fl(~isnan(tbl.fl)));
% plot(yPred);
% plot(tblPredict.fl);
% plot((yPred - min(yPred)) ./ max(yPred - min(yPred)));
% plot((tblPredict.fl - min(tblPredict.fl)) ./ max(tblPredict.fl - min(tblPredict.fl)));
% plot(tblPredict.fwSpeed);
% plot(tblPredict.odorResp);
% legendStr = {'predicted', 'measured'};
yyaxis right
for iWin = 1:numel(odorIntegrationWin)
    varName = ['odorHistory_', num2str(odorIntegrationWin(iWin))];
    plot(tbl.(varName)(~isnan(tbl.fl)));
    legendStr{end + 1} = regexprep(varName, '_', '\\_');
end
% legend(legendStr, 'fontsize', 14)
% legend('predicted', 'measured', 'odorHistory')

% 
%%
% 
% 
% [B, fitInfo] = lasso(table2array(tblShuffle(:, 1:end-1)), tblShuffle.fl, ...
%         'predictorNames', tblShuffle.Properties.VariableNames(1:end-1), 'cv', 10);









