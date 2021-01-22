function regexStr = regexclude(inputStr)
% Generates a regular expression that will fail to find a match (i.e. return empty) when applied to 
% any string that contains the contents of inputStr, but will match any string that does not contain
% inputStr.

% Add escape characters to the input string as needed
specialChars = '\[](){}*+?|^$.<>';
for iChar = 1:numel(specialChars)
    inputStr = strrep(inputStr, specialChars(iChar), ['\', specialChars(iChar)]);
end

% Generate the regex string
if numel(inputStr) > 1
    regexStr = ['^.(?!.*', inputStr, ')(?!', inputStr(2:end), ')|^[^', inputStr(1), ']', ...
        inputStr(2:end), '.*'];
else
    regexStr = ['^[^', inputStr, '](?!.*', inputStr, ')'];
end

end