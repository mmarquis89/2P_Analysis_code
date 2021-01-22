parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\new_PPL201_experiments\Saved_linear_models\'; 
load(fullfile(parentDir, 'nlm_testing_base_standardized.mat'), 'rm');

expNum = 9
histWin =60;

measuredFl = (rm.modelData.fullDataTbl{expNum}.fl);
tbl = rm.modelData.fullDataTbl{expNum};
tblPred = tbl(:, {'moveSpeed', 'odorResp', ['odorHistory_', num2str(histWin)], 'fl'}); 
tblPred{:, 3} = tblPred{:, 3} - min(tblPred{:, 3});
tblFit = tblPred(rm.modelData.fitRowInds{expNum}, :);
tblTest = tblPred(rm.modelData.testRowInds{expNum}, :);

opts = struct();
opts.Display = 'iter';
opts.MaxIter = 2000;
opts.TolFun = 1e-8;
opts.TolX = 1e-8;

modelFun = @nlmFun;
beta0 = ones(1, 6) .* 0.5;
% beta0(5) = -0.5;

mdl = fitnlm(tblFit, modelFun, beta0, 'options', opts)

% Calculate R2 on test dataset
[~, adjR2] = r_squared(tblTest.fl, predict(mdl, tblTest(:, 1:end-1)), mdl.NumCoefficients);

f = figure(1); clf;
f.Color = [1 1 1];
ax(1) = subaxis(3, 1, [1:2], 'ml', 0.05, 'mr', 0.05, 'mb', 0.08, 'sv', 0.08); hold on;
plot(measuredFl);
plot(predict(mdl, tblPred(:, 1:end-1)));
legend({'Measured Fl', 'Predicted Fl'}, 'location' , 'nw');
ax(2) = subaxis(3, 1, 3); hold on;
plot(tblPred{:, 1}, 'linewidth', 1, 'color', 'b');
plot(tblPred{:, 2}, 'linewidth', 1, 'color', 'r');
plot(tblPred{:, 3}, 'linewidth', 3, 'color', 'k');
% odorHistTransformed = (tblPred{:, 3} + 0.00001).^mdl.Coefficients.Estimate(5);
% plot(odorHistTransformed, 'linewidth', 3, 'color', 'm');
linkaxes(ax, 'x');

%% 
cf = mdl.Coefficients.Estimate;
% cf([2 4]) = [1 -1] * 280;
figure(5); clf; hold on;
for i=1:3
%     plot(tblPred{:, i} .* cf(i))
end
plot(tblPred{:, 2} .* (tblPred{:, 3}.^abs(cf(5))));
% plot(cf(4).*tblPred{:,2}.*(tblPred{:, 3} .^ abs(cf(5))))
% plot((cf(4).*tblPred{:,2}.*(tblPred{:, 3} .^ abs(cf(5)))) + (tblPred{:, 2}.*cf(2)))
% plot(measuredFl - cf(6));
%%
a = 0.1;

fcn_1 = @(x) x.^a;
fcn_2 = @(x) 1 - exp(-a.*x);
fcn_3 = @(x) x ./ (1 + x.^a).^(1./a);
% fcn_4 = @(x) log2(x);
% fcn_5 = @(x) log(x);

xx = 0:0.01:3;
figure(1);clf; hold on;
plot(xx, fcn_1(xx));
plot(xx, fcn_2(xx));
plot(xx, fcn_3(xx));
% plot(xx, fcn_4(xx));
% plot(xx, fcn_5(xx));

%%
f = figure(2);clf;
f.Color = [1 1 1];
subaxis(1, 4, 1, 'ml', 0.04, 'mr', 0.04);
histogram(tblPred.moveSpeed, 100);
title('Speed')
subaxis(1, 4, 2);
histogram(tblPred.odorResp, 100);
title('odorResp')
subaxis(1, 4, 3);
histogram(tblPred.(['odorHistory_', num2str(histWin)]), 100);
title('odorHistory')
subaxis(1, 4, 4);
histogram(tblPred.fl, 100);
title('Fl')






