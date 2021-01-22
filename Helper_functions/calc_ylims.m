function ylims = calc_ylims(dataArr, padFactor)
%===================================================================================================
% Calculates and returns y-axis limits suitable for plotting any number in dataArr, encompassing the 
% min and max values as well as some padding on each side. Can optionally specify the size of the 
% padding factor.
%
% INPUTS:
%       dataArr   = numeric array of any dimensions containing the data to be plotted.
%
%       padFactor = <optional> scalar multiple of the absolute range to pad the y-axis limits with.
%                   If an empty array [] is passed the default of 0.1 will be used. 
%
% OUTPUTS:
%     ylims     = numeric vector in the form of [yMin, yMax]. 
%===================================================================================================

if isempty(padFactor)
    padFactor = 0.1;
end

yMax = max(dataArr(:));
yMin = min(dataArr(:));
padding = 0.1 * (yMax - yMin);  % To give a little extra padding
ylims = [yMin - padding, yMax + padding];

end