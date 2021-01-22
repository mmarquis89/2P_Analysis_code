function maxLoc = argmax(inputArr, dim, nanflag)
% Returns the index of the maximum value of a vector or array. If only one argument is provided then
% dim defaults to the lowest dimension with size > 1, and nanflag defaults to 'omitnan'.
if nargin < 3
    nanflag = 'omitnan';
end
if numel(inputArr) == 1
    maxLoc = 1;
else
    if nargin < 2
        sz = size(inputArr);
        dim = find(sz > 1, 1);
    end
    
    [~, maxLoc] = max(inputArr, [], dim, nanflag);
end
end