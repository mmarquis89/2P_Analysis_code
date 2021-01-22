function output = mov_norm(inputVec, winSize)
% Works similarly to movmean() or movmin(), but instead of returning the mean and min, respectively, 
% it normalizes each value to the maximum value in the current window. Only accepts 1-D input 
% vectors.
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
    currWinMax = max(inputVec(currWinStart:currWinEnd));
    output(i) = inputVec(i) / currWinMax;
end

end