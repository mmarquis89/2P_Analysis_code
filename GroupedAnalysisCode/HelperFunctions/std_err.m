function sem = std_err(inputArr, dim, nanflag)
% Calculates standard error of the input array. If only one argument is provided then
% dim defaults to the lowest dimension with size > 1, and nanflag defaults to 'omitnan'.
if nargin < 3
    nanflag = 'omitnan';
end
if nargin < 2
    sz = size(inputArr);
    dim = find(sz > 1, 1);
    if min(sz) < 1 || isempty(dim)
       sem = nan;
       return  
    end
    
    
end

sd = std(inputArr, [], dim, nanflag);
sem = sd ./ (size(inputArr, dim)^0.5);


end