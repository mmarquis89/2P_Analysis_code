 
% ROInames = {'L-SLP', 'L-ANT', 'R-SLP', 'R-ANT'};
ROInames = {'L-SLP', 'L-ANT', 'R-SLP', 'R-ANT', 'Control'};

% ROIDataBaseSub = ROIDataAvg(:, :, 1:4) - repmat(ROIDataAvg(:, :, 5), 1, 1, 4);

stimVols = round(volumeRate * stimEpochs);
% ROIDataBaseSub = ROIDataAvg(stimVols, :, 1:4) - repmat(ROIDataAvg(stimVols, :, 5), 1, 1, 4);
ROIDataBaseSub = ROIDataAvg([1:stimVols(1), stimVols(2):end], :, 1:4) - repmat(ROIDataAvg([1:stimVols(1), stimVols(2):end], :, 5), 1, 1, 4);

ROIDataRs = reshape(ROIDataBaseSub, [], 4);

%%

nBins = 60;

ROISD = std(ROIDataRs);
ROIavg = mean(ROIDataRs);
ROIub = ROIavg + (3.5 * ROISD);

ROIDataTrim = ROIDataRs;





% ROIDataTrim = corrMat;
% ROIDataTrim = lagShiftedCorrMat; 

nROIs = size(ROIDataTrim, 2);
% for i=1:nROIs
%    ROIDataTrim(ROIDataRs(:,i) > ROIub(i), i) = nan; 
% end

[R P RL RU] = corrcoef(ROIDataTrim, 'Rows', 'pairwise');
disp(R)
disp(P)
disp((RU - RL))



f = figure(1);clf;
f.Color = [1 1 1];
for i = 1:nROIs
    for j = 1:nROIs
        
        if i < nROIs && j > 1 && j > i 
        subaxis(nROIs-1, nROIs-1, j-1, i);hold on
        
        histogram2(ROIDataTrim(:,i), ROIDataTrim(:,j), nBins, 'DisplayStyle', 'tile')

        xlabel(ROInames{i});
        ylabel(ROInames{j});
        ax = gca();
        ax.XTickLabel = [];
        ax.YTickLabel = [];
        ax.Title.FontWeight = 'normal';
        ax.Title.Units = 'normalized';
        ax.Title.HorizontalAlignment = 'left';
        ax.Title.Position = [0.05 0.85 0];
        ax.Title.String = ['r = ' num2str(R(i,j), 2)];
        
        suptitle(regexprep(expDate, '\_', '\\_'));
        end
    end
end

