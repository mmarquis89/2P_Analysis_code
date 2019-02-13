
smWin = 3;

ROIDataNorm = [];
for iTrial = 1:nTrials
    ROIDataNorm(:, iTrial, 1) = ROIDataAvg(:,iTrial,1) / max(as_vector(ROIDataAvg(:,iTrial,1)));
end
for iTrial = 1:nTrials
    ROIDataNorm(:, iTrial, 2) = ROIDataAvg(:,iTrial,2) / max(as_vector(ROIDataAvg(:,iTrial,2)));
end

% figure(1);clf;plot(ROIDataNorm(:,:,2)); hold on;
% plot(mean(ROIDataAvg(:,:,1), 2), 'linewidth', 5)
% plot(mean(ROIDataAvg(:,:,2), 2), 'linewidth', 5)

ROIDataLin = as_vector(ROIDataNorm(:,:,1));
ROIDataLin = [ROIDataLin, as_vector(ROIDataNorm(:,:,2))];
figure(1);clf; plot(movmean(ROIDataLin(:,1), smWin));hold on ;plot(movmean(ROIDataLin(:,2), smWin))
xTicks = 1:volumeRate:size(ROIDataLin, 1);
xTickLabels = 1:round(size(ROIDataLin, 1)/volumeRate);
ax = gca; ax.XTick = xTicks; ax.XTickLabel = xTickLabels;