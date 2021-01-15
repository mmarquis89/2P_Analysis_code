function outputData = circ_smooth(inputData, varargin)
% Unwraps the input array and smooths it using the "smoothdata()" function, then re-wraps it and 
% returns the result.
% 
% Inputs: (inputArr, varargin:[range=[0, 2*pi], Dim=first non-singleton, Method='gaussian', 
%                              SmWin, NanFlag='omitnan'])
% Default smooth window is determined heuristically by smoothdata().
% Ignores NaN values in the smoothing window, but will not change any existing NaN values.
%
p = inputParser;
addParameter(p, 'Range', [0 2*pi]);
addParameter(p, 'Dim', []);
addParameter(p, 'Method', 'gaussian');
addParameter(p, 'SmWin', []);
addParameter(p, 'NanFlag', 'omitnan');
parse(p, varargin{:});
range = p.Results.Range;
dim = p.Results.Dim;
method = p.Results.Method;
smWin = p.Results.SmWin;
nanFlag = p.Results.NanFlag;

% Find first non-singleton dimension to operate over if necessary
if isempty(dim)
    while size(inputData, dim) == 1
        dim = dim + 1;
    end
end

% Rescale from 0 to 2*pi
scaledData = (inputData ./ range(2)) * (2*pi);

% Unwrap
uwData = unwrap(scaledData);

% Smooth
if isempty(smWin)
    smUwData = smoothdata(uwData, dim, method, nanFlag);
else
    smUwData = smoothdata(uwData, dim, method, smWin, nanFlag);
end

% Re-wrap
smData = mod(smUwData, 2*pi); 

% Re-scale to original range
outputData = (smData  ./ (2*pi)) .* range(2);

end