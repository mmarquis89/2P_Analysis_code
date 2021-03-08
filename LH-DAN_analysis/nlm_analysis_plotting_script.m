%===================================================================================================
% Load data
%===================================================================================================
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\new_PPL201_experiments\Saved_linear_models\'; 


cellType = 'PPL201';


% Load correct data and set params accordingly
if strcmp(cellType, 'PPL201')
    load(fullfile(parentDir, 'nlm_testing_base_PPL201_normalized_60-60-180.mat'), 'rm');
    load(fullfile(parentDir, 'PPL201_fixedTau_model_results_2.mat'), 'allCoeffs', 'beta0', 'opts');
    histWin = 120;
    all_nOdorStims = [2 1 3 3 4 3 8 3 5]; % Manually selected best one for each experiment
    odorKernelDur = 8;
elseif strcmp(cellType, 'PPL203')
    load(fullfile(parentDir, 'nlm_testing_base_PPL203_normalized_60-30-150.mat'), 'rm');
    load(fullfile(parentDir, 'PPL203_fixedTau_model_results_1.mat'), 'allCoeffs', 'beta0', 'opts');
    histWin = 60;
    all_nOdorStims = 100 * ones(1, 9); % Just use all of them
    odorKernelDur = 18;
end

%% =================================================================================================
% Create source data tables for each experiment
%===================================================================================================
allExpTblPred = {};
try
allVolTimes = {};
for iExp = 1:size(allCoeffs, 2)
    
    expNum = iExp;
    nOdorStims = all_nOdorStims(iExp);
    
    %===============================================================================================
    % Set up model variables
    %===============================================================================================
    histTerm = ['odorHistory_', num2str(histWin)];
    tbl = rm.modelData.fullDataTbl{expNum};
    tblPred = tbl(:, {'moveSpeed', 'odorResp', histTerm, 'fl'});
    if expNum == 8 && odorKernelDur == 8
        tblPred.fl(5245:end) = nan;
    end
    
    % Find first n odor onset times
    volTimes = rm.sourceData.volTimes{expNum};
    odorVols = rm.sourceData.odorVols{expNum};
    odorOnsetVols = regexp(regexprep(num2str(odorVols), ' ', ''), '(?<=0)1');
    avgWinStartVols = odorOnsetVols;
    volumeRate = rm.sourceData.volumeRate(expNum);
    kernelLenVols = floor(volumeRate * odorKernelDur);
    avgWinEndVols = avgWinStartVols + kernelLenVols - 1;
    warning('off');
    
    % Calculate new odor response kernel
    odorRespKernel = zeros(kernelLenVols, 1);
    avgCount = 0;
    stimCount = 1;
    while stimCount < numel(avgWinStartVols) && avgCount < nOdorStims
        if avgWinEndVols(stimCount) < avgWinStartVols(stimCount + 1)
            currKernelInds = avgWinStartVols(stimCount):(avgWinEndVols(stimCount));
            if ~any(isnan(tblPred.fl(currKernelInds)))
                odorRespKernel = odorRespKernel + tblPred.fl(currKernelInds);
                avgCount = avgCount + 1;
            end
        end
        stimCount = stimCount + 1;
    end
    odorRespKernel = odorRespKernel ./ avgCount;
    odorRespKernel = odorRespKernel - odorRespKernel(1);
    originalKernel = rm.sourceData.odorRespFilter{expNum};
    if sum(odorRespKernel) == 0
        odorRespKernel = originalKernel;
    end
    
    % Convolve with odor onset vols to get new predicted response to each odor
    odorOnsetVector = zeros(size(odorVols));
    odorOnsetVector(odorOnsetVols) = 1;
    odorRespVector = conv(odorOnsetVector, odorRespKernel);
    odorRespVector = odorRespVector(1:numel(odorOnsetVector));
    
    % Add a steady drift term and replace odorResp term;
    tblPred.volsFromExpStart = normalize(1:numel(volTimes))';
    tblPred.odorResp = odorRespVector';
    tblPred = tblPred(:, [1 2 3 5 4]);
    
    % Fit and subtract the drift term
    modelSpec = ['fl ~ volsFromExpStart'];
    mdl = fitlm(tblPred, modelSpec);
    predFl = predict(mdl, tblPred(:, 1:end-1));
    tblPred.fl = tblPred.fl - predFl;
    
    % Create final source data table
    tblPred.odorOnsets = odorOnsetVector';
    tblPred = tblPred(:, [2 1 6 3 5]);
    
    allExpTblPred{expNum} = tblPred;
    allVolTimes{iExp} = volTimes;
end
catch ME; rethrow(ME); end


%% =================================================================================================
% Fit new models without an odor history term
%===================================================================================================
noOdorHistMdls = {}; allR2_noOdorHist = [];
for iExp = 1:numel(allExpTblPred)
    currTbl = allExpTblPred{iExp};
    
    beta0_noOdorHist = beta0([1:4, 6]);
    tbl_noOdorHist = currTbl(:, [1:3, 5]);
    modelFun = @nlmFun_no_odorHist;
    
    noOdorHistMdls{iExp} = fitnlm(tbl_noOdorHist(~isnan(tbl_noOdorHist.fl), :), modelFun, ...
            beta0_noOdorHist, 'options', opts);
    allR2_noOdorHist(iExp) = noOdorHistMdls{iExp}.Rsquared.Adjusted;
end

%% =================================================================================================
% Fit new models without any adaptation terms
%===================================================================================================
noAdaptationMdls = {}; allR2_noAdaptation = [];
for iExp = 1:numel(allExpTblPred)
    currTbl = allExpTblPred{iExp};
    
    beta0_noAdaptation = beta0([3:6]);
    tbl_noAdaptation = currTbl(:, [1:2, 4:5]);
    modelFun = @nlmFun_no_adaptation;
    
    noAdaptationMdls{iExp} = fitnlm(tbl_noAdaptation(~isnan(tbl_noAdaptation.fl), :), modelFun, ...
            beta0_noAdaptation, 'options', opts);
    allR2_noAdaptation(iExp) = noAdaptationMdls{iExp}.Rsquared.Adjusted;
end

%% =================================================================================================
% Fit new models without any adaptation OR odor hist terms
%===================================================================================================
minimalMdls = {}; allR2_minimal = [];
for iExp = 1:numel(allExpTblPred)
    currTbl = allExpTblPred{iExp};
    
    beta0_minimal = beta0([3 4 6]);
    tbl_minimal = currTbl(:, [1 2 5]);
    modelFun = @nlmFun_minimal;
    
    minimalMdls{iExp} = fitnlm(tbl_minimal(~isnan(tbl_minimal.fl), :), modelFun, ...
            beta0_minimal, 'options', opts);
    allR2_minimal(iExp) = minimalMdls{iExp}.Rsquared.Adjusted;
end
%% =================================================================================================
% Fit new models with only an odor term
%===================================================================================================
odorOnlyMdls = {}; allR2_odorOnly = [];
for iExp = 1:numel(allExpTblPred)
    currTbl = allExpTblPred{iExp};
    
    beta0_odorOnly = beta0([1 6]);
    tbl_odorOnly = currTbl(:, [1 5]);
    modelFun = @nlmFun_odorOnly;
    
    odorOnlyMdls{iExp} = fitnlm(tbl_odorOnly(~isnan(tbl_odorOnly.fl), :), modelFun, ...
            beta0_odorOnly, 'options', opts);
    allR2_odorOnly(iExp) = odorOnlyMdls{iExp}.Rsquared.Adjusted;
end
%% =================================================================================================
% Fit new models with only an moveSpeed term
%===================================================================================================
speedOnlyMdls = {}; allR2_speedOnly = [];
for iExp = 1:numel(allExpTblPred)
    currTbl = allExpTblPred{iExp};
    
    beta0_speedOnly = beta0([2 6]);
    tbl_speedOnly = currTbl(:, [2 5]);
    modelFun = @nlmFun_speedOnly;
    
    speedOnlyMdls{iExp} = fitnlm(tbl_speedOnly(~isnan(tbl_speedOnly.fl), :), modelFun, ...
            beta0_speedOnly, 'options', opts);
    allR2_speedOnly(iExp) = speedOnlyMdls{iExp}.Rsquared.Adjusted;
end

%% Plot comparison of R2 values across different model versions

figDir = parentDir;
saveFig = 1;
fileNameSuffix = '_noOdorHist-noAdapt';
% fileNameSuffix = '_noAdapt-noOdorHist';

allR2 = allCoeffs{'adjR2', :};

% Exclude experiments with models that failed to fit
if strcmp(cellType, 'PPL203')
    failedFits = [2 3];
    
%     r2comp = [allR2_minimal; allR2_noAdaptation; allR2_noOdorHist; allR2];
%     labels = {'base', '-adapt', '-odorHist', 'full'};
%     

    labels = {'base', '-odorHist', '-adapt', 'full'};
    r2comp = [allR2_minimal; allR2_noOdorHist; allR2_noAdaptation; allR2];
    
elseif strcmp(cellType, 'PPL201')
    failedFits = 1;
%     
%     r2comp = [allR2_minimal; allR2_noAdaptation; allR2_noOdorHist; allR2];
%     labels = {'base', '-adapt', '-odorHist', 'full'};
    

    labels = {'base', '-odorHist', '-adapt', 'full'};
    r2comp = [allR2_minimal; allR2_noOdorHist; allR2_noAdaptation; allR2];
%     
end
r2comp(:, failedFits) = [];

% Plot figure
f = figure(3); clf;
f.Color = [1 1 1];
f.Position(3:4) = [450 530];
plot(r2comp, '-o', 'color', [0.5 0.5 0.5]);
hold on; 
errorbar(1:size(r2comp, 1), mean(r2comp, 2), std_err(r2comp, 2), 's', 'color', 'k', 'linewidth', 2);
xlim([0.5, size(r2comp, 1) + 0.5]);
ax = gca();
ax.XTick = 1:size(r2comp, 1);
ax.XTickLabel = labels;
title(cellType);
ylabel('Adj R2');
ax.FontSize = 12;

if saveFig
    [p, ~, stats] = anova1(r2comp', [], 'off');
    c = multcompare(stats, 'ctype', 'bonferroni', 'display', 'off');
    c = c(:, [1 2 6]);
    f.UserData.pVals = c;
    fileName = ['nlm_R2_comparison_', cellType, fileNameSuffix];
    save_figure(f, figDir, fileName);
end


%% Plot example data from a selected experiment

expNum = 9;

timeWin = [250 750];
% timeWin = [];

figDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\new_PPL201_experiments\Saved_linear_models\Thesis figs';
saveFig = 1;
fileNameSuffix = '';

% Get data for current experiment
currTbl = allExpTblPred{expNum};
currCoeffs = allCoeffs{1:6, expNum};
[predFl, mdlOdor, mdlMoveSpeed] = nlmFun_double_tau_fixed(currCoeffs, currTbl{:, 1:end-1});
rawMoveSpeed = currTbl.moveSpeed;
rawOdor = currTbl.odorResp;
odorHist = currTbl.(histTerm);
odorOnsets = currTbl.odorOnsets;
measFl = currTbl.fl;
volTimes = allVolTimes{expNum};

f = figure(3); clf;
f.Color = [1 1 1];
f.Position(3:4) = [1900 600];

clear ax;
ax(1) = subaxis(2, 1, 1)
hold on;
plot(volTimes, measFl, 'color', rgb('orange'));
plot(volTimes, predFl, 'color', rgb('green'));
ax(1).FontSize = 14;
expID = allCoeffs.Properties.VariableNames{expNum}; 
title(expID);
ax(1).YTickLabel = [];
if ~isempty(timeWin)
    xlim(timeWin);
end

ax(2) = subaxis(2, 1, 2);
hold on;
plot(volTimes, odorHist, 'color', 'k');
plot(volTimes, mdlOdor, 'color', 'r');
plot(volTimes, rawMoveSpeed, 'color', 'm');
plot(volTimes, mdlMoveSpeed, 'color', 'b');
plot(volTimes, odorOnsets, 'color', 'g');
plot(volTimes, rawOdor, 'color', 'c');
ax(2).FontSize = 14;
ax(2).YTickLabel = [];
if ~isempty(timeWin)
    xlim(timeWin);
end
linkaxes(ax, 'x');

if saveFig
    fileName = ['nlm_example_data_', cellType, '_', expID, '_', num2str(timeWin(1)), '-', ...
            num2str(timeWin(2)), fileNameSuffix];
    save_figure(f, figDir, fileName);
end



%%

expNum = 2;

coeffs = allCoeffs{:, expNum}';
[yhat, ~, ~] = nlmFun_double_tau_fixed(coeffs, allExpTblPred{expNum}{:, 1:end-1});
volTimes = allVolTimes{expNum};

f = figure(100);clf;
f.Color = [1 1 1];
ax(1) = subaxis(1, 1, 1);
hold on;
plot(volTimes, allExpTblPred{expNum}.fl);
plot(volTimes, yhat)
legend({'Measured', 'Predicted'}, 'location', 'nw', 'fontsize', 12)
    xlim([0 volTimes(end)])
ylabel('Coeffs\_1', 'FontSize', 16);

% coeffs = allCoeffs_1{:, expNum}';
% [yhat, ~, ~] = nlmFun_double_tau_fixed(coeffs, temp{:, 1:end-1});

% ax(2) = subaxis(2, 1, 2);   
% hold on;
% plot(volTimes, temp.fl);
% plot(volTimes, yhat)
% legend({'Measured', 'Predicted'}, 'location', 'nw', 'fontsize', 12)
%     xlim([0 volTimes(end)])
% ylabel('Coeffs\_2', 'FontSize', 16);
% 
% linkaxes(ax);

