function outputArr = repeat_smooth(inputArr, nReps, varargin)
% Repeatedly applies the "smoothdata()" function to an input array
% Inputs: (inputArr, nReps, varargin:[Dim=1, Method='gaussian', SmWin])
% Default smooth window is determined heuristically by the smoothdata().

p = inputParser;
addParameter(p, 'Dim', 1);
addParameter(p, 'Method', 'gaussian');
addParameter(p, 'SmWin', []);
parse(p, varargin{:});
dim = p.Results.Dim;
method = p.Results.Method;
smWin = p.Results.smWin;

outputArr = inputArr;

for iRep = 1:nReps
    if isempty(smWin)
        outputArr = smoothdata(outputArr, dim, method);
    else
        outputArr = smoothdata(outputArr, dim, method, smWin);
    end
end
end