function overwrite = confirm_save(saveFile, varargin)
% Checks to see if a file exists and if it does, asks user whether they are sure they want to 
% overwrite it.

% Parse optional argument
p = inputParser;
addParameter(p, 'DialogText', 'This will overwrite an existing file...continue?');
parse(p, varargin{:});
dlgText = p.Results.DialogText;

% Check if file exists and confirm overwrite
overwrite = 1;
if exist(saveFile, 'file')
    dlgAns = questdlg(dlgText, 'Warning', 'Yes', 'No', 'No');
    if strcmp(dlgAns, 'No')
        overwrite = 0;
    end
end

end