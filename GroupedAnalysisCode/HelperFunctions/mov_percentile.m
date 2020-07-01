function output = mov_percentile(inputVec, winSize, percentile)
% Works similarly to movmean() or movmin(), but instead of returning the mean and min, respectively, 
% returns an arbitrary percentile of that window. Only accepts 1-D input vectors.
output = [];
halfWin = round(winSize / 2);
for i = 1:numel(inputVec)
    currWinStart = i - halfWin;
    currWinEnd = i + halfWin;
    if currWinStart < 1
        currWinStart = 1;
    end
    if currWinEnd > numel(inputVec)
       currWinEnd = numel(inputVec); 
    end
    currWinVals = sort(inputVec(currWinStart:currWinEnd));
    pctIndex = round(numel(currWinVals) .* percentile);
    if pctIndex < 1
        pctIndex = 1;
    end
    output(i) = currWinVals(pctIndex);
end

end