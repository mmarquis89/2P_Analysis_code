


currTrial = 16;


currVectAvgWedge = bl.dffVectAvgWedge(:, currTrial);
barPos = bl.panelsPosX;



vol2pFrame = sample_lookup(bl.panelsDisplayRate, bl.volumeRate);
volFrames = vol2pFrame.convert(1:numel(bl.volTimes));
volFrames(volFrames > numel(bl.panelsFrameTimes)) = numel(bl.panelsFrameTimes);

figure(1); clf; hold on
plot(bl.volTimes, bl.dffVectAvgRad(:, currTrial));
plot(bl.volTimes, (barPos(volFrames) * ((2*pi)/100)) - pi);

figure(2);clf;hold on
xPositions = circshift(-180:3.75:(180 - 3.75), -5);
plotX = barPos(volFrames);
newPlotX = [];
for i=1:numel(xPositions)
    newPlotX(plotX == i) = xPositions(i);
end


plot(newPlotX, bl.dffVectAvgRad(:, currTrial), '.')