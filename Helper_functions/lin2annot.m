function outputArr = lin2annot(linArr, outputDims)
%===================================================================================================
% 
% Simple function for converting a vector of concatenated trial data into an intuitively-formatted
% 2D array. The output array will be formatted such that each row represents one trial and each 
% column represents one time point in that trial. Reading such an array from left-to-right and 
% top-to-bottom (as if reading written text) would result in the input to this function (linArr).
% 
% See "annot2lin" for more info.
% 
% INPUTS: 
%       linArr     = the input array of linear annotation/event data, formatted as described above
%       outputDims = the size of the linear array's 2D counterpart with dims: [nTrials, nTimePoints] 
%
% OUTPUTS:
%       outputArr  = the 2D output array with dimensions as described above
%
%===================================================================================================

outputArr = reshape(linArr, fliplr(outputDims))';

end