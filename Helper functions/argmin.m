function minLoc = argmin(inputArr, dim, nanflag)
% Returns the index of the minimum value of a vector or array. If only one argument is provided then
% dim defaults to the lowest dimension with size > 1, and nanflag defaults to 'omitnan'.
if nargin < 3
    nanflag = 'omitnan';
end
if nargin < 2
    sz = size(inputArr);
    dim = find(sz > 1, 1);
end
[~, minLoc] = min(inputArr, [], dim, nanflag);

end