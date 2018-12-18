function outputArr = repeat_smooth(inputArr, smWin, dim, nReps)
% Repeatedly applies the "movmean()" function to an input array

outputArr = inputArr;
for iRep = 1:nReps
    outputArr = movmean(outputArr, smWin, dim, 'omitnan'); 
end
end