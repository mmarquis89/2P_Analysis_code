function dffData = calc_dFF(respData, baselineData, avgDims)
%================================================================================================
% Calculates averaged dF/F from baseline and response period data. Can average across either 
% trials, volumes, or both. Always averages across both when calculating the baseline. Any 'inf' 
% values that are caused by dividing by zero will be replaced with zeros in the output array.
%
% INPUTS:
%       respData     = data from the period you want to get dF/F for: [y, x, plane, volume, trial]
%
%       baselineData = data in the same format, but from the baseline period only
%
%       avgDims      = a vector containing the dimensions you want to average over. Valid values 
%                      include either [5] for trial-averaging, [4] for volume averaging, or [4 5] 
%                      to average across both trials and volumes.
%
% OUTPUTS:
%       dffData = trial-averaged dF/F data, averaged across avgDims
%
%================================================================================================

avgDims = sort(avgDims, 'descend');

respMean = respData;
for iDim = 1:length(avgDims)
    respMean = squeeze(mean(respMean, avgDims(iDim)));                                  % --> [y, x, plane, volume] (ex. for avgDims = 5)
end
baselineMean = squeeze(mean(mean(baselineData, 5), 4));                                 % --> [y, x, plane]
baselineMeanRep = repmat(baselineMean, 1, 1, 1, size(respMean, 4), size(respMean, 5));  % --> [y, x, plane, volume]

dffData = (respMean - baselineMeanRep) ./ baselineMeanRep;                              % --> [y, x, plane, volume]

dffData(isinf(dffData)) = 0;
end