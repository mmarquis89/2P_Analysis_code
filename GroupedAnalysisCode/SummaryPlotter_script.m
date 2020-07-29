expList = load_expList();
expList = expList(contains(expList.expName, 'D-ANT'), :);
expList = expList(~contains(expList.groupName, 'newExpts'), :);

plt = SummaryPlotter2D(expList);

%% Save SummaryPlotter object
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Saved_summaryPlotter2D_objects';
fileName = 'SummaryPlotter.mat'; 
save(fullfile(parentDir, fileName), 'plt');

%%
maxVals = [1500 1500 1500 1500];
minVals = [0 0 0 0];

expNum = 61;

plotVarNames = {'TypeD', 'fwSpeed'};
maxVals = [2000 200 2000, 200, 200, 200];
minVals = [-1000 0 0 0 0 0];
spacerRowOverride = 0;
ignoreSavedParams = 0;

saveParams = 0;
saveFig = 0;
fileNameTag = 'D-ANT';

SV = 0.06;
SH = 0.045;
ML = 0.04;
MR = 0.02;
MT = 0.05;
MB = 0.1;

disp(expList.expID{expNum});
nPlots = numel(plotVarNames);
if nPlots ~= 4
    spDims = numSubplots(nPlots);
else
    spDims = [2 2];
end
plt.params.ignoreSavedParams = ignoreSavedParams;
f = figure(1);
f.Color = [1 1 1];
for iPlot = 1:nPlots
    ax = subaxis(spDims(1), spDims(2), iPlot, ...
            'sv', SV, 'sh', SH, 'ml', ML, 'mr', MR, 'mt', MT, 'mb', MB);
    plt.params.expID = expList.expID{expNum};
    plt.params.plotVarName = plotVarNames{iPlot};
    
    % Set min and max values for current plot data
    if ~isempty(maxVals)
        plt.params.maxVal = maxVals(iPlot);
    end
    if ~isempty(minVals)
        plt.params.minVal = minVals(iPlot);
    end
    
    % Make plot and save colormap range if desired
    plt.plot_summary(ax);
    if saveParams
        plt = plt.save_colormap_range();
    end
    
    % Remove redundant internal axis labels
    if spDims(1) > 1 && iPlot <= (nPlots - spDims(2))
        ax.XLabel.String = '';
    end
    if spDims(2) > 1 && ~ismember(iPlot, 1:(spDims(2)):nPlots)
        ax.YLabel.String = '';
    end
    
    % Add expDate to title of first plot
    if iPlot == 1
       ax.Title.String = [expList.expID{expNum}, ' — ', ax.Title.String]; 
    end
    
    % Adjust colormap
    cm = viridis([], 0.02);
    if strcmp(expList.groupName{expNum}, 'gaplessAcq') && ~spacerRowOverride
        cm = [cm(1:end-1, :); 0 0 0];
    end
    colormap(cm);
    
    ax.XTick = ax.XTick(1:2:end);
    ax.XTickLabel = ax.XTickLabel(1:2:end, :);
    
end

% Save figure
if saveFig
    saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Figs';
    fileName = ['2D_summary_plots_', fileNameTag, '_', expList.expID{expNum}];
    save_figure(f, saveDir, fileName)
end


