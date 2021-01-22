function outputStr = capitalize(str)
% Capitalizes just the first character of the input string. Returns nan if passed an empty string.
if numel(str) > 1
    outputStr = [upper(str(1)), str(2:end)];
elseif numel(str) == 1
    outputStr = upper(str);
else
    outputStr = nan;
end
end