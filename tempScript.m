%%

figure; hold on
plot(1:10, medTimes);
plot(1:10, smallTimes);
plot(1:10, tinyTimes);


medPx = repmat(131072, 1, 10);
smallPx = repmat(32768, 1, 10);
tinyPx = repmat(8192, 1, 10);


x1 = repmat([1:10]', 3, 1);
x2 = [repmat(131072, 10, 1); repmat(32768, 10, 1); repmat(8192, 10, 1)];
y = [medTimes'; smallTimes'; tinyTimes'];

figure(1); clf; hold on
models = []; coeffs = [];
for iSize = 1:size(testTimes, 1)
    currTimes = testTimes(iSize, :)';
    currModel = fit([1000:1000:5000]', currTimes, 'poly2');
    models{iSize} = currModel;
    coeffs(iSize, :) = [currModel.p1, currModel.p2, currModel.p3];
    fitData = [0:100:20000];
    plotData = currModel.p1 .* (fitData.^2) + currModel.p2 .* fitData + currModel.p3;
    plot([1000:1000:5000], currTimes, 'o');
    plot(fitData, plotData);
end










