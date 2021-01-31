parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\new_PPL201_experiments\Saved_linear_models\'; 
% load(fullfile(parentDir, 'nlm_testing_base_PPL201_standardized_120-60-240.mat'), 'rm');
% load(fullfile(parentDir, 'nlm_testing_base_PPL203_standardized_60-30-240.mat'), 'rm');
load(fullfile(parentDir, 'nlm_testing_base_PPL201_normalized_60-60-180.mat'), 'rm');
% load(fullfile(parentDir, 'nlm_testing_base_PPL203_normalized_60-30-150.mat'), 'rm');

% Choose experiment and odor history window size
expNum = 9;   % 1-9
histWin = 120; % 120, 180, 240

all_nOdorStims = [2 1 3 3 4 3 8 3 5]; % Manually selected best one for each experiment
nOdorStims = all_nOdorStims(expNum);
% nOdorStim = 100; % Calculates average of first n odor stims as the response kernal

odorKernelDur = 8;

%===================================================================================================
% Set up model variables
%===================================================================================================
histTerm = ['odorHistory_', num2str(histWin)];
tbl = rm.modelData.fullDataTbl{expNum};
tblPred = tbl(:, {'moveSpeed', 'odorResp', histTerm, 'fl'}); 
if expNum == 8
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

% Plot new kernel overlaid with the orignal one and the mean movement speed
f = figure(1);clf; 
f.Color = [1 1 1];
hold on;
originalKernel = rm.sourceData.odorRespFilter{expNum};
if sum(odorRespKernel) == 0
    odorRespKernel = originalKernel;
end
meanOdorSpeed = rm.sourceData.odorRespMeanSpeed{expNum};
plot(volTimes(1:numel(odorRespKernel)), odorRespKernel/max(abs(odorRespKernel), [], 'omitnan'), 'linewidth', 2);
plot(volTimes(1:numel(originalKernel)), originalKernel/max(abs(originalKernel), [], 'omitnan'), 'linewidth', 2);
xL = xlim();
plot((1:numel(meanOdorSpeed)) ./ xL(2), (meanOdorSpeed./max(abs(meanOdorSpeed), [], 'omitnan')), 'linewidth', 2);
plot(xL, [0 0], '--', 'color' ,'k', 'linewidth', 2)
xlim(xL);
legend({'New odorResp kernel', 'Original odorResp kernel', 'Avg moveSpeed (original kernel)'}, ...
        'fontsize', 10, 'location', 'best')

% Convolve with odor onset vols to get new predicted response to each odor
odorOnsetVector = zeros(size(odorVols));
odorOnsetVector(odorOnsetVols) = 1;
odorRespVector = conv(odorOnsetVector, odorRespKernel);
odorRespVector = odorRespVector(1:numel(odorOnsetVector));

% Add a steady drift term and replace odorResp term;
tblPred_2 = tblPred;
tblPred_2.volsFromExpStart = normalize(1:numel(volTimes))';
tblPred_2.odorResp = odorRespVector';
tblPred_2 = tblPred_2(:, [1 2 3 5 4]);

% Fit and subtract the drift term
modelSpec = ['fl ~ volsFromExpStart'];
mdl = fitlm(tblPred_2, modelSpec);
predFl = predict(mdl, tblPred_2(:, 1:end-1));
tblPred_2.fl = tblPred_2.fl - predFl;

%===================================================================================================
% Fit a nonlinear model to remaining data and plot results
%===================================================================================================
tblPred_3 = tblPred_2(:, {'odorResp', 'fl'});
tblPred_3.odorOnsets = odorOnsetVector';
tblPred_3 = tblPred_3(:, [1 3 2]);
opts = struct();
opts.Display = 'iter';
opts.MaxIter = 200;
opts.DerivStep = .001;
opts.TolFun = 1e-8;
opts.TolX = 1e-8;

temp = tblPred_3;
temp.moveSpeed = tblPred_2.moveSpeed;
temp.odorHistory_120 = tblPred_2.odorHistory_120;
% temp.odorHistory_60 = tblPred_2.odorHistory_60;
temp = temp(:, [1 4 2 5 3]);
% 
% modelFun = @nlmFun_double_tau_fit;
% beta0 = [0 0 3000 500 1.5 0.5 -0.1 0]; % [f_odor, f_speed, tau_odor, tau_speed, b_odor, b_speed, b_odorHist, intercept]

modelFun = @nlmFun_double_tau_fixed;
beta0 = [0.1 0.1 1 0.5 -0.1 0]; % [f_odor, f_speed, b_odor, b_speed, b_odorHist, intercept]

mdl = fitnlm(temp(~isnan(temp.fl), :), modelFun, beta0, 'options', opts);

coefTbl = mdl.Coefficients(:,1);
% coefTbl.Row = {'f_odor', 'f_moveSpeed', 'tau_odor', 'tau_moveSpeed', 'b_odor', 'b_moveSpeed', 'b_odorHistory', 'intercept'};
coefTbl.Row = {'f_odor', 'f_moveSpeed', 'b_odor', 'b_moveSpeed', 'b_odorHistory', 'intercept'};
if ~exist('allCoeffs', 'var')
    allCoeffs = [coefTbl; table(mdl.Rsquared.Adjusted, 'variableNames', {'Estimate'}, 'rownames', {'adjR2'})];
    allCoeffs.Properties.VariableNames = (rm.modelData.expID(expNum));
end
allCoeffs.(rm.modelData.expID{expNum}) = [coefTbl.Estimate; mdl.Rsquared.Adjusted];
disp(coefTbl);

% [yhat, yhat_odor, yhat_moveSpeed] = nlmFun_double_tau_fit(mdl.Coefficients.Estimate', ...
%         temp{:, 1:end-1});
[yhat, yhat_odor, yhat_moveSpeed] = nlmFun_double_tau_fixed(mdl.Coefficients.Estimate', ...
        temp{:, 1:end-1});

if ~all(isnan(yhat_odor))
    % Plot final results
    clear ax;
    nPlots = 3;
    f = figure(3);clf;
    f.Color = [1 1 1];
    ax(1) = subaxis(nPlots, 1, 1, 'ml', 0.04, 'mr', 0.04, 'mb', 0.04, 'mt', 0.06, 'sv', 0.03);
    hold on
    plot(volTimes, temp.fl);
    plot(volTimes, yhat);
    legend({'Measured', 'Predicted'}, 'location', 'nw', 'fontsize', 12)
    xlim([0 volTimes(end)])
    t = title([rm.modelData.expID{expNum}, '  —  Final model adjR2: ', ...
            num2str(mdl.Rsquared.Adjusted, 2)]);
    t.FontSize = 12;
    ax(2) = subaxis(nPlots, 1, 2);
    hold on;
    plot(volTimes, yhat_odor, 'color', 'r', 'linewidth', 1)
    plot(volTimes, yhat_moveSpeed, 'color', 'b');
    plot(volTimes, temp.(histTerm), 'color', 'k', 'linewidth', 3);
    xlim([0 volTimes(end)])
    legend({'moveSpeed', 'odorResp', regexprep(histTerm, '\_', '\\_')}, ...
            'location', 'nw', 'fontsize', 12);
    ax(3) = subaxis(nPlots, 1, 3);
    hold on;
    plot(volTimes, temp.fl, 'linewidth', 1);
    plot(volTimes, temp.fl - yhat, 'linewidth', 1, 'color', 'k');
    plot([volTimes(1), volTimes(end)], [0 0], '--', 'linewidth', 2, 'color', 'r')
    legend({'Measured Fl', 'Unexplained by model'}, 'location', 'nw', 'fontsize', 12)
    
    xlim([0 volTimes(end)])
    linkaxes(ax);
    
end


