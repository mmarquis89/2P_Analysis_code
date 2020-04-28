function outputArr = multi_mean(inputArr, dims)
% ==================================================================================================
% Computes the mean over multiple dimensions of an array and returns the squeezed result
% ==================================================================================================
outputArr = inputArr; 
dims = sort(dims); % Ensures that dimensions will stay in the same places during averaging
for iDim = 1:numel(dims)
    outputArr = mean(outputArr, dims(iDim));
end 
outputArr = squeeze(outputArr);

end