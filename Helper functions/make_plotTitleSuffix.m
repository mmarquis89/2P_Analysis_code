function plotTitleSuffix = make_plotTitleSuffix(groupNames)
plotTitleSuffix = ['  —  ', groupNames{1}];
for iStim = 1:numel(groupNames)
    if iStim > 1
        plotTitleSuffix = [plotTitleSuffix, ' vs. ', groupNames{iStim}];
    end
end
plotTitleSuffix = [plotTitleSuffix, ' (top to bottom)'];
end