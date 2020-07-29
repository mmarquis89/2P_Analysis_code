
smWin = 50;
nReps = 40;

aData = a.analysisOutput;
slidingWinSize = aData.paramSnapshot.slidingWinSize;

base = mov_percentile(aData.rawFlData, slidingWinSize, 0.05)';
flData = (aData.rawFlData - base) ./ base;


rsFl = reshape(flData, [], 180);
figure(18); clf;
imagesc(rsFl')
colormap(viridis)


f = figure(19);clf; 
f.Color = [1 1 1];
subaxis(3, 1, 1); hold on; 
plot(aData.rawFlData(~isnan(aData.rawFlData)), 'color', 'b')
plot(repeat_smooth(aData.rawFlData(~isnan(aData.rawFlData)), 40, 'SmWin', 50), ...
        'linewidth', 3, 'color', 'r')
title('Raw Fl')

subaxis(3, 1, 2); hold on; 
plot(flData(~isnan(flData)), 'color', 'b')
plot(repeat_smooth(flData(~isnan(flData)), 40, 'SmWin', 50), 'linewidth', 3, 'color', 'r')
title('Sliding dF/F')

subaxis(3, 1, 3); hold on; 
plot(aData.fwSpeedData(~isnan(aData.fwSpeedData)), 'color', 'b')
plot(repeat_smooth(aData.fwSpeedData(~isnan(aData.fwSpeedData)), 40, 'SmWin', 50), 'linewidth', 3, 'color' , 'r')
title('fw speed')


