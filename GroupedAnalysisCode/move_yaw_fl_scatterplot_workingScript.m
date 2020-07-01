testSpeed = a.analysisOutput.moveSpeedData;
testYaw = a.analysisOutput.yawSpeedData;
testFl = a.analysisOutput.rawFlData;

intFl = round(testFl) - min(round(testFl)) + 1;
cm = jet(round(max(intFl)));
cm = cm(intFl, :) * 0.9;
%
figure(15);clf; scatter(testSpeed, testYaw, 10, cm, 'o')
xlabel('speed');
ylabel('yaw speed');


% Plot 2D tiled speed vs. yaw plot colored by fluorescence 

nBins = 30;
figure(16);clf;
h = histogram2(testSpeed, testYaw, nBins, 'displaystyle', 'tile');
xEdges = h.XBinEdges;
yEdges = h.YBinEdges;

binMeans = [];
for x = 1:numel(xEdges)-1
    for y = 1:numel(yEdges)-1
        xMin = xEdges(x);
        xMax = xEdges(x + 1);
        yMin = yEdges(y);
        yMax = yEdges(y + 1);
        binMeans(x, y) = mean(testFl(testSpeed >= xMin & testSpeed < xMax & ...
                               testYaw >=yMin & testYaw < yMax), 'omitnan');
%         binMeans(x, y) = mean(testFl();
    end
end
imagesc(flipud(binMeans'));
set(gca, 'XTick', 1:numel(xEdges));
set(gca, 'XTickLabel', xEdges);
set(gca, 'YTick', 1:numel(yEdges));
set(gca, 'YTickLabel', yEdges(end:-1:1));



