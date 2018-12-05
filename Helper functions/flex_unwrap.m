function outputData = flex_unwrap(inputData, varargin)
%===================================================================================================
% 
% Works essentially the same as the built-in "unwrap()" function, except that it can also be used 
% on datasets that wrap at arbitrary endpoints (instead of just 0 and 2*pi like the original). The
% default functionality should be essentially identical to the original function. 
% 
% INPUTS:
%       inputData = the time series data that you want to unwrap 
% 
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       'Range' = (default: [min(inputData(:)), max(inputData(:)]) the values at which the data 
%                  wraps, in the form of a two element vector [min max]
%
%       'Tolerance' = (default: midpoint of inputData range) the minimum jump magnitude to be 
%                     treated as a wrap
%
%       'Dimension' = (default: 1) if array has > 1 non-singleton dimension, the dimension that it 
%                     should be unwrapped along
%
%===================================================================================================

% Parse optional arguments
p = inputParser;
addParameter(p, 'Range', [min(inputData(:)), max(inputData(:))]);
addParameter(p, 'Tolerance', 0.5 * abs(max(inputData(:)) - min(inputData(:))));
addParameter(p, 'Dimension', 1);
parse(p, varargin{:});
range = p.Results.Range;
tolerance = p.Results.Tolerance;
uwDim = p.Results.Dimension;

wrapMag = max(range) - min(range); % get the total size of the wrapping range
nDims = sum(size(inputData) > 1); % get number of non-singleton dimensions of inputData

if nDims == 1
    outputData = local_unwrap(inputData, wrapMag, tolerance);    
else
    % Bring target dimension to the 1st position
    permOrder = [uwDim:nDims, 1:uwDim - 1];
    permData = permute(inputData, permOrder);
    
    % Reshape input data into a matrix
    sz = size(permData);
    resData = reshape(permData, [sz(1), prod(sz(2:end))]);
    
    % Unwrap each column along target dimension
    uwResData = zeros(size(resData));
    for iCol = 1:size(resData, 2)
        uwResData(:,iCol) = local_unwrap(resData(:,iCol), wrapMag, tolerance);
    end
    
    % Return output to original shape
    outputData = reshape(uwResData, sz);
    outputData = ipermute(outputData, permOrder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Local Function  %%%%%%%%%%%%%%%%%%%
    function localOut = local_unwrap(resData, wrapMag, tolerance)
        
        wrapCount = 0;
        
        for iSamp = 1:(size(resData) - 1)
            
            % Bring sample up to current unwrapped value
            localOut(iSamp,:) = resData(iSamp) + (wrapMag * wrapCount);
            
            if((abs(resData(iSamp + 1) - resData(iSamp))) > (abs(tolerance)))  %if diff is greater than tolerance, increment or decrement wrapCount
                
                if resData(iSamp + 1) < resData(iSamp)   % if the phase jump is negative, increment wrapCount
                    wrapCount = wrapCount + 1;
                else             % if the phase jump is positive, decrement wrapCount
                    wrapCount = wrapCount - 1;
                end
            end
        end
        localOut(iSamp + 1,:) = resData(iSamp + 1) + (wrapMag * wrapCount); % add wrapMag * wrapCount to last element
    end

end