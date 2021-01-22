function tbl = hd(A, k)
% Wrapper for the built-in head() function that changes the default number of rows that are 
% displayed from 8 to 1
if nargin < 2
    k = 1;
end
tbl = head(A, k);
end