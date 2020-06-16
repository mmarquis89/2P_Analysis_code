function [binMidpoints, binMeans, binSEM] = bin_data(xData, yData, nBins, minX)
        
        plotDataTable = table(xData, yData, 'variableNames', {'speed', 'fl'});
        plotDataSort = sortrows(plotDataTable, 'speed');
        plotDataSort = plotDataSort(plotDataSort.speed > minX, :);
        binSize = round(size(plotDataSort, 1) / nBins);
        binEdges = 1:binSize:(size(plotDataSort, 1) - 1);
        flBinMeans = []; SEM = []; speedBinMidpoints = [];
        for iBin = 1:numel(binEdges)
            if iBin < numel(binEdges)
                binInds = binEdges(iBin):(binEdges(iBin + 1) - 1);
            else
                binInds = binEdges(iBin):size(plotDataSort, 1);
            end
            flBinMeans(iBin) = mean(plotDataSort{binInds, 'fl'}, 'omitnan');
            SEM(iBin) = std_err(plotDataSort{binInds, 'fl'});
            speedBinMidpoints(iBin) = plotDataSort{binInds(round(numel(binInds) / 2)), 'speed'};
        end
        binMidpoints = speedBinMidpoints;
        binMeans = flBinMeans;
        binSEM = SEM;
end