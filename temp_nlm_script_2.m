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
% nOdorStims = 100; % Calculates average of first n odor stims as the response kernal

odorKernelDur = 8;

%===================================================================================================
% Set up model variables
%===================================================================================================
histTerm = ['odorHistory_', num2str(histWin)];
tbl = rm.modelData.fullDataTbl{expNum};
tblPred = tbl(:, {'moveSpeed', 'odorResp', histTerm, 'fl'}); 

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
measuredFl = tblPred_2.fl;


%===================================================================================================
% Fit best model possible without using the odor response and plot results
%===================================================================================================
% Fig 2 subplot descriptions (top to bottom):
% 1) Measured fluorescence overlaid with initial (sans odor-response data) linear model prediction
% 2) Overlaid predictor variables for the current experiment 
% 3) Measured - predicted (AKA blue line in subplot #1 - orange line in subplot #1)
% 4) Nonlinear model's predicted odor response overlaid onto data from previous subplot
% 5) The estimated odor response vector before and after running it through the adaptation model
% 6) The remainting unexplained fluorescence after subtracting the final model prediction

% modelSpec = ['fl ~ moveSpeed + ', histTerm, ' + moveSpeed:volsFromExpStart'];
modelSpec = ['fl ~ moveSpeed * ', histTerm]; 
mdl = fitlm(tblPred_2, modelSpec)
predFl = predict(mdl, tblPred_2(:, 1:end-1));

% Plot results
clear ax;
nPlots = 6;
f = figure(2);clf;
f.Color = [1 1 1];
ax(1) = subaxis(nPlots, 1, 1, 'ml', 0.04, 'mr', 0.04, 'mb', 0.04, 'mt', 0.06, 'sv', 0.03); 
hold on
plot(volTimes, measuredFl);
plot(volTimes, predFl);
legend({'measured', 'predicted'}, 'location', 'nw', 'fontsize', 10)
xlim([0 volTimes(end)])
ax(2) = subaxis(nPlots, 1, 2);
hold on;
plot(volTimes, tblPred_2.moveSpeed);
plot(volTimes, tblPred_2.odorResp);
plot(volTimes, tblPred_2.(histTerm));
xlim([0 volTimes(end)])
legend({'moveSpeed', 'odorResp', regexprep(histTerm, '\_', '\\_')}, ...
        'location', 'nw', 'fontsize', 10);

% Subtract all non-odor resp predictions out of the model and plot remaining data
ax(3) = subaxis(nPlots, 1, 3);
hold on;
plot(volTimes, measuredFl - predFl);
xlim([0 volTimes(end)])
legend('unexplained #1', 'location', 'nw', 'fontsize', 10)
linkaxes(ax, 'x');

%===================================================================================================
% Fit a nonlinear model to remaining data and plot results
%===================================================================================================
tblPred_3 = tblPred_2(:, {'odorResp', 'fl'});
tblPred_3.odorOnsets = odorOnsetVector';
tblPred_3 = tblPred_3(:, [1 3 2]);
tblPred_3.fl = tblPred_3.fl - predFl;
opts = struct();
opts.Display = 'iter';
opts.MaxIter = 100;
opts.TolFun = 1e-6;
opts.TolX = 1e-6;

modelFun = @nlmFun_tau_fixed;
beta0 = [0.8 0 0]; % Starting coeff values (fixed tau)

% modelFun = @nlmFun_tau_fit;
% beta0 = [0.8 5000 0 0]; % Starting coeff values (fit tau)

mdl = fitnlm(tblPred_3(~isnan(tblPred_3.fl), :), modelFun, beta0, 'options', opts)
% tau(expNum) = mdl.Coefficients{'b2', 'Estimate'};

% Plot odorResp modeling results
ax(4) = subaxis(nPlots, 1, 4); hold on
plot(volTimes, measuredFl - predFl);
plot(volTimes, predict(mdl, tblPred_3(:, 1:end-1)));
legend({'unexplained #1', 'predicted'}, 'location', 'nw', 'fontsize', 10);
xlim([0 volTimes(end)])

ax(5) = subaxis(nPlots, 1, 5);
hold on;
plot(volTimes, tblPred_3.odorResp);
plot(volTimes, predict(mdl, tblPred_3(:, 1:end-1)));
xlim([0 volTimes(end)]);
legend({'Estimated odor response', 'Adapted odor response'}, 'location', 'nw', 'fontsize', 10);

ax(6) = subaxis(nPlots, 1, 6);
hold on;
plot(volTimes, (measuredFl - predFl) - predict(mdl, tblPred_3(:, 1:end-1)));
% plot(volTimes, predict(mdl, tblPred_3(:, 1:end-1)));
xlim([0 volTimes(end)]);
legend('Unexplained (final)', 'location', 'nw', 'fontsize', 10)
linkaxes(ax, 'x');

%===================================================================================================
% Use the result as a new estimated odor response and re-fit the linear model to the full dataset
%===================================================================================================
% Fig 3 subplot descriptions (top to bottom):
% 1) Measured fluorescence overlaid with final model predicted fluorescence
% 2) Overlaid predictor variables for the final model
% 3) Measured fluorescence overlaid with the remaining unexplained data

tblPred_4 = tblPred_2;
tblPred_4.odorResp = predict(mdl, tblPred_3(:, 1:end-1));
% modelSpec = ['fl ~ moveSpeed + odorResp +', histTerm, '+ moveSpeed:volsFromExpStart']; 
modelSpec = ['fl ~ moveSpeed + odorResp + ', histTerm, ' + moveSpeed:', histTerm]; 
mdl = fitlm(tblPred_4, modelSpec)

% Plot final results
clear ax;
nPlots = 3;
f = figure(3);clf;
f.Color = [1 1 1];
ax(1) = subaxis(nPlots, 1, 1, 'ml', 0.04, 'mr', 0.04, 'mb', 0.04, 'mt', 0.06, 'sv', 0.03); 
hold on
plot(volTimes, measuredFl);
plot(volTimes, predict(mdl, tblPred_4(:, 1:end-1)));
legend({'Measured', 'Predicted'}, 'location', 'nw', 'fontsize', 12)
xlim([0 volTimes(end)])
t = title([rm.modelData.expID{expNum}, '  —  Final model adjR2: ', ...
        num2str(mdl.Rsquared.Adjusted, 2)]);
t.FontSize = 12;
ax(2) = subaxis(nPlots, 1, 2);
hold on;
plot(volTimes, tblPred_4.moveSpeed, 'color', 'b');
plot(volTimes, tblPred_4.odorResp, 'color', 'r', 'linewidth', 1);
plot(volTimes, tblPred_4.(histTerm), 'color', 'k', 'linewidth', 3);
xlim([0 volTimes(end)])
legend({'moveSpeed', 'odorResp', regexprep(histTerm, '\_', '\\_')}, ...
        'location', 'nw', 'fontsize', 12);
ax(3) = subaxis(nPlots, 1, 3);
hold on;
plot(volTimes, measuredFl, 'linewidth', 1);
plot(volTimes, measuredFl - (predict(mdl, tblPred_4(:, 1:end-1))), 'linewidth', 1, 'color', 'k');
plot([volTimes(1), volTimes(end)], [0 0], '--', 'linewidth', 2, 'color', 'r')
legend({'Measured Fl', 'Unexplained by model'}, 'location', 'nw', 'fontsize', 12)

xlim([0 volTimes(end)])
linkaxes(ax);

% adjR2_fit(expNum) = mdl.Rsquared.Adjusted;
% adjR2_fixed_3000(expNum) = mdl.Rsquared.adjusted;

disp(['Tau = ' num2str(tau(expNum), 5)]);
disp(num2str(tau, 5));




