function write_to_log(msg, idStr)

if nargin < 2
   idStr = '';
else
    idStr = [' (', idStr, ')'];
end

try
    % Appends lines to a log file for overnight data processing
    myFile = fopen('/home/mjm60/dataProcessingLog.txt', 'a');
    fprintf(myFile, [datestr(datetime), idStr, ':  ', msg '\r\n']);
    fclose(myFile);
catch
end