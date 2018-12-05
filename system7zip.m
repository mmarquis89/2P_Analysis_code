function system7zip(parentDir, archiveName, archiveType, filterString, test)
%=====================================================================================================================
% ADD FILES TO A COMPRESSED ARCHIVE
% Calls 7zip from the Windows command line to archive files in one of several formats. Archive file will be located
% in the parent directory. After completion the function optionally attempts to test the integrity of the archive, 
% again using a system command to run 7zip, and saves the results of the test as a .txt file in the archive directory. 
%
% Inputs:
%       parentDir    = the directory containing the files you want to archive (e.g. 'D:\2p Behavior Video\2017_07_30')
%       archiveName  = the name (minus file extension) of the archive file to be created (e.g. 'RawFrames')
%       archiveType  = string containing the name of any 7zip-supported archive format (e.g. 'zip')
%                       The valid values for archiveType are: 'zip', '7z', 'tar', and 'wim'
%       filterString = string to specify which files in parentDir are to be added to the archive (e.g. *sid_0_tid*_mp4)
%       test         = boolean indicating whether to test the integrity of the archive after creating it
%======================================================================================================================

% Make sure archiveType is valid
validTypes = {'zip', '7z', 'tar', 'wim'};
assert(sum(strcmp(archiveType, validTypes)) == 1, ['Error: "', archiveType, ...
        '" is not a supported archive type. The valid types are: "zip", "7z", "tar", and "wim"']);

% Make sure a file with this name does not already exist
fullFileName = [archiveName, '.', archiveType];
assert(exist(fullfile(parentDir, fullFileName), 'file')== 0, 'Error: an archive with that name already exists in this directory');

% Add files to archive
disp(['Adding files to ', fullFileName, '...'])
cmdStr_1 = [parentDir(1:2), ' & ']; % Switch working directory to the correct drive if necessary
cmdStr_2 = ['cd ', parentDir, ' & ']; % Set working directory to parentDir
cmdStr_3 = '"C:\Program Files\7-Zip\7z.exe"'; % Need to provide the entire path to 7zip.exe for some reason
cmdStr_4 = [' a ', fullFileName, ' ', filterString, ' & ']; % Create 7zip command
cmdStr_5 = [' t ', fullFileName, ' >', fullFileName, '_log.txt & ']; % Archive testing command
if test
    system([cmdStr_1, cmdStr_2, cmdStr_3, cmdStr_4, cmdStr_3, cmdStr_5]);
else
    system([cmdStr_1, cmdStr_2, cmdStr_3, cmdStr_4]);
end

end%function