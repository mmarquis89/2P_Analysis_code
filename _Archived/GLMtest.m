

binSizeSec = 4;

binSizeVols = round(binSizeSec * mD(1).volumeRate);
dffData = mD(4).expDffMat(:, 1);
Y = dffData(binSizeVols:end, 1);

% Get stim data vector
stimData = zeros(size(dffData));
stimTiming = mD(1).optoStimTiming;
odorOnsetTimes = stimTiming(1):sum(stimTiming(2:3)):mD(1).trialDuration;
odorOffsetTimes = odorOnsetTimes + stimTiming(2);
if odorOffsetTimes(end) == mD(1).trialDuration
    odorOnsetTimes(end) = [];
    odorOffsetTimes(end) = [];
end 
sec2vol = sample_lookup(mD(1).volumeRate, 1);
odorOnsetVols = sec2vol.convert(odorOnsetTimes);
odorOffsetVols = sec2vol.convert(odorOffsetTimes);
for iStim = 1:numel(odorOnsetVols)
    stimData(odorOnsetVols(iStim):odorOffsetVols(iStim)) = 1;
end

% Create design matrix
X = [];
for iVol = 1:(numel(stimData) - binSizeVols + 1)
    X(iVol, :) = stimData(iVol:iVol + binSizeVols - 1);
end


b = glmfit(X, smoothdata(Y, 'gaussian', 5), 'normal');