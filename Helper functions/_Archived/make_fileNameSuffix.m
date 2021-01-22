function fileNameSuffix = make_fileNameSuffix(stimNames)
fileNameSuffix = ['_', stimNames{1}];
for iStim = 1:numel(stimNames)
    if iStim > 1
        fileNameSuffix = [fileNameSuffix, 'vs', stimNames{iStim}];
    end
end
end