function varargout = split_vector(inputVec)
%===================================================================================================
% 
% Returns each element of a numeric input vector as an individual output argument.
%
% INPUTS:
%       inputVec = numeric vector that you want to split into individual variables.
%
% OUTPUTS:
%       varargout = each element of the input vector is returned as a separate variable.
%
%===================================================================================================

varargout = num2cell(inputVec);

end