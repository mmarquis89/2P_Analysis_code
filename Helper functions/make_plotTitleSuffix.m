function plotTitleSuffix = make_plotTitleSuffix(stimNames)
plotTitleSuffix = ['  —  ', stimNames{1}];
for iStim = 1:numel(stimNames)
    if iStim > 1
        plotTitleSuffix = [plotTitleSuffix, ' vs. ', stimNames{iStim}];
    end
end
plotTitleSuffix = [plotTitleSuffix, ' (top to bottom)'];
end