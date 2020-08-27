function x = exp_floor(n, startingExp)
% Recursively rounds a number down to the next integer power of 10
% Starts at 10e10 by default, but can pass a second argument to make it work with larger numbers
if nargin < 2
    startingExp = 10;
end
if n < 10^startingExp
    x = exp_floor(n, startingExp - 1);
else
    x = 10^startingExp;
end
end