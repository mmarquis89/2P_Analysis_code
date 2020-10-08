function save_figure(figHandle, saveDir, fileName, overwriteFlag)
%===================================================================================================
% Wrapper function for saving figures as a .png .pdf, and .fig file that asks to make sure the user
% wants to overwrite a file of the same name if one exists. The .fig file is saved within a 
% separate subdirectory (which is created if necessary) called "figFiles".
%
% INPUTS:
%
%       figHandle     = handle to the figure to be saved
%
%       saveDir       = directory where you want the .png files saved
%
%       fileName      = output file name
%
%       overwriteFlag = (optional positional) pass "skipConfirmation" to skip overwrite confirmation
%
%===================================================================================================

% Create save directory if necessary
if ~isdir(saveDir)
    mkdir(saveDir);
end

% Warn user and offer to cancel save if this will overwrite existing files
if nargin > 3 && strcmp(overwriteFlag, 'skipConfirmation')
    overwrite = 1;
else
    overwrite = 1;
    if exist(fullfile(saveDir, [fileName, '.fig']), 'file') ~= 0 || ...
            exist(fullfile(saveDir, [fileName, '.png']), 'file') ~= 0
        dlgAns = questdlg(['Saving this figure will overwrite one or more existing files in this ', ...
            'directory...are you sure you want to do this?'], 'Warning', 'Yes', 'No', 'No');
        if strcmp(dlgAns, 'No')
            overwrite = 0;
            disp('Saving cancelled')
        end
    end
end

% Save figure files
if overwrite 
    export_fig(fullfile(saveDir, fileName), '-png', '-nocrop', figHandle);
    s = warning('error', 'export_fig:transparency');
    try
        export_fig(fullfile(saveDir, fileName), '-pdf', '-nocrop', figHandle);
    catch
        disp('Transparency detected...using print() instead of export_fig() to save PDF')
        print(figHandle, '-dpdf', '-bestfit', fullfile(saveDir, fileName));
    end
    warning(s);
    if ~isdir(fullfile(saveDir, 'figFiles'))
        mkdir(fullfile(saveDir, 'figFiles'))
    end
    savefig(figHandle, fullfile(saveDir, 'figFiles', fileName));
end

end