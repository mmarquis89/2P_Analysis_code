
goodFlyNums = [3     7     8    10    15    17    18    19];
analysisData = allData;
goodFlyData = analysisData(ismember([analysisData.flyNum], goodFlyNums));
for i = 1:numel(goodFlyData)
    goodFlyData(i).flyNum = find(goodFlyNums == goodFlyData(i).flyNum);
end
nFlies = numel(goodFlyNums);
analysisData = goodFlyData;

analysisData = allData;
badFlyNums = 1:24; 
badFlyNums = badFlyNums(~ismember(badFlyNums, goodFlyNums));
badFlyData = analysisData(ismember([analysisData.flyNum], badFlyNums));
for i = 1:numel(badFlyData)
    badFlyData(i).flyNum = find(badFlyNums == badFlyData(i).flyNum);
end
nFlies = numel(badFlyNums);
analysisData = badFlyData;

runningFlyNums = [1 2 11 13 20 24];
analysisData = allData;
runningFlyData = analysisData(ismember([analysisData.flyNum], runningFlyNums));
for i = 1:numel(runningFlyData)
    runningFlyData(i).flyNum = find(runningFlyNums == runningFlyData(i).flyNum);
end
nFlies = numel(runningFlyNums);
analysisData = runningFlyData;