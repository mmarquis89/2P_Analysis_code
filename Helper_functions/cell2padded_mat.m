% Takes a cell array of variable-length column vectors and converts them to a matrix by padding all
% vectors with NaN to match the size of the longest one
function outputMat = cell2padded_mat(inputCell)
    maxLen = max(cellfun(@numel, inputCell));
    outputMat = [];
    for i = 1:numel(inputCell)
        if numel(inputCell{i}) == maxLen
            outputMat(:, i) = inputCell{i};
        else
            padVec = nan(maxLen - numel(inputCell{i}), 1);
            outputMat(:, i) = [inputCell{i}; padVec];
        end        
    end
end