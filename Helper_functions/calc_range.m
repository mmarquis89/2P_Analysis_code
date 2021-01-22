function range = calc_range(dataArr, varargin)
%==============================================================================================================
% Calculates the maximum absolute value of an array and returns [min, max] values for the array that
% are symmetrically centered at zero. Useful for example in calculating ranges for plots using a 
% colormap centered at zero. Can optionally multiply output by a scalar value or use the standard 
% deviation in place of the max for the calculation.
%
% Note: optional arguments can be supplied either positionally or as name/value pairs, with the function
%       assuming that any scalar is 'scalingFactor' and any string is 'method'.
%
% INPUTS:
%       dataArr        = numeric array of any dimension containing the data to be encompassed by the range.
%
%       scalingFactor  = <OPTIONAL> scalar value to multiply the output range by. Default is 1.
%
%       method         = <OPTIONAL> string specifying which metric to use to calculate the range. 
%                        Default is "max" but can be changed to "stdDev".
%
% OUTPUTS:
%       range   = numeric vector in the form of [min, max] containing the range values (e.g., [-3, 3]) 
%===============================================================================================================

% Parse and validate arguments
p = inputParser;
addOptional(p, 'optPos_1', []);
addOptional(p, 'optPos_2', [], @ischar);
addParameter(p, 'scalingFactor', 1);
addParameter(p, 'calcMethod', 'max', @ischar);
parse(p, varargin{:});
scalingFactor = p.Results.scalingFactor;
calcMethod = p.Results.calcMethod;

% If optional positional arguments were provided, figure out which is which
if ~isempty(p.Results.optPos_1)
   optArg = p.Results.optPos_1;
   if isscalar(optArg)
       scalingFactor = optArg;
   elseif ischar(optArg)
       calcMethod = optArg;
   end
end
if ~isempty(p.Results.optPos_2)
   optArg = p.Results.optPos_2;
   if isscalar(optArg)
       scalingFactor = optArg;
   elseif ischar(optArg)
       calcMethod = optArg;
   end
end

if strcmpi(calcMethod, 'max')
    
    % Use max value for range calculation
    minMax = [];
    minMax(1) = min(dataArr(:));
    minMax(2) = max(dataArr(:));
    range(2) = max([abs(min(minMax)), max(minMax)]);
    range(1) = -range(2);
    range = range * scalingFactor;

elseif strcmpi(calcMethod, 'stdDev')
    
   % Use standard deviation for range calculation
    stdDev = std(dataArr(:));
    range = [-stdDev, stdDev] .* scalingFactor;
    
end
end