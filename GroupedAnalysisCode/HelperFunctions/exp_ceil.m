function exp_ceil(n, startingExp)
% Recursively rounds a number up to the next integer power of 10
% Starts at 10e-5 by default, but can pass a second argument to make it work with smaller numbers
if nargin < 2
    startingExp = -5;
end
if n > 10^startingExp
    x = exp_ceil(n, startingExp + 1);
else
    x = 10^startingExp;
end
end