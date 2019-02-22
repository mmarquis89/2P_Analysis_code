 
% ROInames = {'L-SLP', 'L-ANT', 'R-SLP', 'R-ANT'};
% ROInames = {'L-SLP', 'L-ANT', 'R-SLP', 'R-ANT', 'Control'};

flData = ROIDffAvg;

allROIs = size(flData, 3);
currROIs = allROIs - 1;
stimVols = round(volumeRate * stimEpochs);

ROIData = flData(:, :, 1:currROIs);
baselineData = repmat(flData(:, :, end), 1, 1, currROIs);
ROIDataBaseSub = ROIData;%- baselineData;

% ROIDataBaseSub = ROIDataBaseSub(stimVols, :, :);
ROIDataBaseSub = ROIDataBaseSub([1:stimVols(1), stimVols(2):end], :, :);

ROIDataRs = reshape(ROIDataBaseSub, [], size(ROIDataBaseSub, 3));

%%


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

