function outputArr = repeat_smooth(inputArr, nReps, varargin)
% Repeatedly applies the "smoothdata()" function to an input array
% Inputs: (inputArr, nReps, varargin:[Dim=1, Method='gaussian', SmWin, NanFlag='omitnan'])
% Default smooth window is determined heuristically by the smoothdata().
% Ignores NaN values in the smoothing window, but will not change any existing NaN values.

p = inputParser;
addParameter(p, 'Dim', 1);
addParameter(p, 'Method', 'gaussian');
addParameter(p, 'SmWin', []);
addParameter(p, 'NanFlag', 'omitnan');
parse(p, varargin{:});
dim = p.Results.Dim;
method = p.Results.Method;
smWin = p.Results.SmWin;
nanFlag = p.Results.NanFlag;

outputArr = inputArr;

for iRep = 1:nReps
    if isempty(smWin)
        outputArr = smoothdata(outputArr, dim, method, nanFlag);
    else
        outputArr = smoothdata(outputArr, dim, method, smWin, nanFlag);
    end
end

outputArr(isnan(inputArr)) = nan;

end