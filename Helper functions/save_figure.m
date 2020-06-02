function save_figure(figHandle, saveDir, fileName)
%===================================================================================================
% Wrapper function for saving figures as a .png .pdf, and .fig file that asks to make sure the user
% wants to overwrite a file of the same name if one exists. The .fig file is saved within a 
% separate subdirectory (which is created if necessary) called "figFiles".
%
% INPUTS:
%
%       figHandle   = handle to the figure to be saved
%
%       saveDir     = directory where you want the .png files saved
%
%       fileName    = output file name
%
%===================================================================================================

% Create save directory if necessary
if ~isdir(saveDir)
    mkdir(saveDir);
end

% Warn user and offer to cancel save if this will overwrite existing files
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

% Save figure files
if overwrite
    export_fig(fullfile(saveDir, fileName), '-png', '-pdf', '-nocrop', figHandle);
    if ~isdir(fullfile(saveDir, 'figFiles'))
        mkdir(fullfile(saveDir, 'figFiles'))
    end
    savefig(figHandle, fullfile(saveDir, 'figFiles', fileName));
end

end