function sem = std_err(inputArr, dim, nanflag)
% Calculates standard error of the input array. Defaults to dim=1 and nanflag='omitnan' if only one 
% argument is provided.

if nargin < 3
    nanflag = 'omitnan';
end
if nargin < 2
   dim = 2; 
end

sd = std(inputArr, [], dim, nanflag);
sem = sd ./ (size(plotData, dim)^0.5);


end