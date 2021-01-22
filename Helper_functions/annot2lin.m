function outputArr = annot2lin(inputArr)
%===================================================================================================
% 
% Simple function for quickly linearizing 2D annotation or trial data in an intuitive way. The input 
% array should be formatted such that each row represents one trial and each column represents one 
% time point in that trial. Reading such an array from left-to-right and top-to-bottom (as if 
% reading written text) would result in the output vector, which simply consists of all time points 
% throughout the experiment in chronological order.
% 
% Use "lin2annot" to reverse this process.
% 
% INPUTS: 
%       inputArr = the input array of annotation/event data, formatted as described above
%
% OUTPUTS:
%       outputArr = the linearized version of the input array
%
%===================================================================================================

outputArr = inputArr';
outputArr = outputArr(:)';

end