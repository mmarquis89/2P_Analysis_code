

% Reorder file/folder names with timestamp first so I can sort them in the order they were created 

parentDir = 'B:\Dropbox (HMS)\2P Data\Imaging Data\2017_10_22';

% B:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_06
%


allFolders = dir(fullfile(parentDir, '*sid*tid*'));

folderNames = {allFolders.name}';

newNames = [];
for iFolder = 1:length(folderNames)
   fName = folderNames{iFolder};
   sidIdx = strfind(fName, '_sid');
   sidTid = fName(sidIdx+1:end-4);
   timeStamp = fName(sidIdx - 15:sidIdx-1);
   stimInfo = fName(7:sidIdx - 17);
   newNames{iFolder} = [timeStamp, '_', sidTid, '_', stimInfo, '.tif'];
end

count = 1;
for iFolder = 1:length(folderNames)
    disp(num2str(count));
    oldFolder = fullfile(parentDir, folderNames{iFolder});
    newFolder = fullfile(parentDir, newNames{iFolder});
    movefile(oldFolder, newFolder,'f');
    count = count + 1;
end
