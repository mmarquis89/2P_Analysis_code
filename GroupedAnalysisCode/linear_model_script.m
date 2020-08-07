parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';

% Set parameters
expID = '20190315-3';
roiName = 'TypeD-R';
trialNums = [1];
maxSpeed = 100;
smWinVols = 3;
smWinFrames = 3;
smReps = 10;
odorIntegrationWin = [30];
ftLagVols = 3;

% -------------------- Load data -------------------------% 
% Load exp and trial metadata
[expMd, trialMd] = load_metadata({expID}, parentDir);
trialMd = trialMd(trialNum,:);
md = innerjoin(expMd, trialMd);
volTimes = calc_volTimes(md.nVolumes, md.volumeRate, md.trialDuration, md.originalTrialCount);

% Load and smooth imaging rawFl data
roiData = load_roi_data({expID}, parentDir);
roiData = roiData(strcmp(roiData.roiName, roiName), :);
roiData = roiData(trialNum, :);
fl = smoothdata(roiData.rawFl{:}, 'gaussian', smWinVols);

% Load and smooth FicTrac data, then downsample to match volume rate
FRAME_RATE = 25;
ft = load_ft_data({expID}, parentDir);
ft = ft(trialNum, :);
fwSpeed = repeat_smooth(ft.fwSpeed{:} .* 4.5 .* FRAME_RATE, smReps, 'smWin', smWinFrames);
fwSpeed(fwSpeed > maxSpeed) = maxSpeed;
yawSpeed = abs(repeat_smooth(rad2deg(ft.yawSpeed{:}) .* FRAME_RATE, smReps, 'smWin', smWinFrames));
fwSpeedVols = [];
yawSpeedVols = [];
for iVol = 1:numel(volTimes)
    dsFrame = argmin(abs(ft.frameTimes{:} - volTimes(iVol)));
    fwSpeedVols(iVol) = fwSpeed(dsFrame);
    yawSpeedVols(iVol) = yawSpeed(dsFrame);
end

% Apply lag to FicTrac data 
fwSpeedVols = fwSpeedVols([ones(1, ftLagVols), 1:end-ftLagVols]);
yawSpeedVols = yawSpeedVols([ones(1, ftLagVols), 1:end-ftLagVols]);

% Load odor event data and create odor event vector
eventData = load_event_data({expID}, parentDir);
odorEventData = eventData.odor;
odorCommandVols = odorEventData.create_logical_array(md.nVolumes, volTimes, table({expID}, 1, ...
        'variableNames', {'expID', 'trialNum'}));
    
% Calculate integrated odor values

odorIntegrationWinVols = round(odorIntegrationWin * md.volumeRate);
integratedOdor = [];
for iWin = 1:numel(odorIntegrationWinVols)
    currIntegrationWin = odorIntegrationWinVols(iWin);
    for iVol = 1:numel(odorCommandVols)
        if iVol <= currIntegrationWin
            integratedOdor(iWin, iVol) = sum(odorCommandVols(1:iVol));
        else
            integratedOdor(iWin, iVol) = sum(odorCommandVols(iVol - currIntegrationWin:iVol));
        end
    end
end
% ------------------- Convolution ------------------------%

% Create odor response filter
gaussian = @(a, b, c, x) a .* exp(-( (x - b).^2 ./ (2 * c^2) ));
% a1 = 0.3;
% b1 = 1 .* md.volumeRate;
% c1 = 0.3 .* md.volumeRate;
% a2 = 1;
% b2 = 3 .* md.volumeRate;
% c2 = 1 .* md.volumeRate;
a1 = 0.2;
b1 = 0.6 .* md.volumeRate;
c1 = 0.2 .* md.volumeRate;
a2 = 1;
b2 = 2.5 .* md.volumeRate;
c2 = 0.8 .* md.volumeRate;
g1 = gaussian(a1, b1, c1, 1:1:(round(6 * md.volumeRate)));
g2 = gaussian(a2, b2, c2, 1:1:(round(6 * md.volumeRate)));

figure(1);clf; hold on;
plot(volTimes(1:(round(6 * md.volumeRate))), g1);
plot(volTimes(1:(round(6 * md.volumeRate))), g2);
plot(volTimes(1:(round(6 * md.volumeRate))), g1 + g2, 'linewidth', 3);
plot_stim_shading([0,4])
xlim([-2, 4])
odorRespFilter = [0, g1 + g2];

% Extract onset volumes from odor command vector
odorOnsetVols = (regexp(regexprep(num2str(odorCommandVols'), ' ', ''), '01')) + 1;
odorOnsetVector = zeros(size(odorCommandVols));
odorOnsetVector(odorOnsetVols) = 1;

% Convolve with odor response filter to get predicted fl responses
odorRespVector = conv(odorOnsetVector, odorRespFilter);
odorRespVector = odorRespVector(1:numel(odorOnsetVector));

% Generate model

% Create table of predictor and response variables
volumesFromTrialStart = 1:numel(fl);
tbl = table(fwSpeedVols', yawSpeedVols', odorRespVector, volumesFromTrialStart', ...
        'VariableNames', {'fwSpeed', 'yawSpeed', 'odorResp', 'volsFromStart'});
for iWin = 1:size(integratedOdor, 1) 
    varName = ['odorHistory_', num2str(odorIntegrationWin(iWin))];
    tbl.(varName) = integratedOdor(iWin, :)';
end    
tbl.fl = fl;

if ~isnan(md.pmtShutoffVols{:})
    tbl.fl(md.pmtShutoffVols{:}) = nan;
end

% Scale all variables so that the max value is 1;    
for iCol = 1:size(tbl, 2)
    currData = tbl{:, iCol};
    tbl{:, iCol} = tbl{:, iCol} ./ max(tbl{:, iCol});
%     tbl{:, iCol} = tbl{:, iCol} - mean(tbl{:, iCol}, 'omitnan');% Center all variables at zero
end




%     
% figure(2); clf; hold on; plot(sort(tbl.fwSpeed)); plot(sort(tbl.yawSpeed)); ... 
%         plot(sort(tbl.odorResp)); plot(sort(tbl.odorHistory)); plot(sort(tbl.fl))
% legend(tbl.Properties.VariableNames)

% Split data into 100 trials and add every other trial to test vs. predict datasets
splitInd = 1:size(tbl, 1);
rsSplitInd = reshape(splitInd, md.nVolumes / md.originalTrialCount, []);
fitInds = as_vector(rsSplitInd(:, 1:2:end));
predictInds = as_vector(rsSplitInd(:, 2:2:end));
tblFit = tbl(fitInds, :);
tblPredict = tbl(predictInds, :);
% tblFit = tbl(predictInds, :);
% tblPredict = tbl(fitInds, :);


% Drop rows with invalid fluorescence values
tblFit = tblFit(~isnan(tblFit.fl), :);
tblPredict = tblPredict(~isnan(tblPredict.fl), :);

% tblShuffle = tbl(randperm(size(tbl, 1)), :);
% midpoint = floor(size(tblShuffle, 1) / 2);
% tblFit = tblShuffle(1:midpoint, :);
% tblPredict = tblShuffle((midpoint + 1):end, :);


% Fit model
mdl = fitlm(tblFit(:, [1:3, 5:6]))
varNames = tbl.Properties.VariableNames([1:3, 5:6]);
for i=1:(numel(varNames) - 1)
    figure(i + 16); clf; hold on; plot(tbl.(varNames{i}), tbl.fl, '.');
    p = polyfit(tbl.(varNames{i})(~isnan(tbl.fl)), tbl.fl(~isnan(tbl.fl)), 1);
    f = polyval(p, tbl.(varNames{i}), 'omitnan'); plot(tbl.(varNames{i}), f, '-');
    legend('Raw data', [num2str(round(p(1), 2)), 'x + ', num2str(round(p(2), 2))])
    xlabel(varNames{i}); ylabel('fl'); title('Raw data scatter')
    figure(i + 1);clf; plotAdded(mdl, varNames{i}); 
    figure(i + 12);clf; plotAdjustedResponse(mdl, varNames{i});
end
figure(1);clf; plotEffects(mdl);

% figure; plotAdded(mdl, 'odorHistory_30');

% 
% 
% 
% [yPred, yCI] = predict(mdl, tblPredict(:, 1:end-1));
% figure(2);clf;hold on; 
% plot(yPred);
% plot(tblPredict.fl);
% % plot((yPred - min(yPred)) ./ max(yPred - min(yPred)));
% % plot((tblPredict.fl - min(tblPredict.fl)) ./ max(tblPredict.fl - min(tblPredict.fl)));
% % plot(tblPredict.fwSpeed);
% % plot(tblPredict.odorResp);
% legendStr = {'predicted', 'measured'};
% for iWin = 1:numel(odorIntegrationWin)
%     varName = ['odorHistory_', num2str(odorIntegrationWin(iWin))];
%     plot(tblPredict.(varName));
%     legendStr{end + 1} = regexprep(varName, '_', '\\_');
% end
% legend(legendStr, 'fontsize', 14)
% % legend('predicted', 'measured', 'odorHistory')
% 
% 
% 





