
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\allExperiments';
saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';

% Get list of all experiments in parent directory
expMdFiles = dir(fullfile(parentDir, '*_expMetadata.csv'));
fileNames = {expMdFiles.name};
expList = regexp(fileNames, '.*(?=_expMetadata\.csv)', 'match', 'once')';

% Find list of unique files after removing expID prefix
allFiles = dir(parentDir);
allFileNames = {allFiles.name}';
allFileNames = allFileNames(cellfun(@(x) numel(x) > 2, allFileNames));
fileNameSuffixes = regexp(allFileNames, '(?<=\d{8}-\d_).*', 'match', 'once');
uniqueFiles = unique(fileNameSuffixes(~cellfun(@isempty, fileNameSuffixes)));

% Initialize checklist table
fileChecklistArr = zeros(numel(expList), numel(uniqueFiles));
fileChecklistTable = array2table(fileChecklistArr, 'VariableNames', regexprep(uniqueFiles, ...
        '\.', '_'));

% Fill in checklist table
for iExp = 1:numel(expList)
    for iFile = 1:numel(uniqueFiles)
        currFileName = [expList{iExp}, '_', uniqueFiles{iFile}];
        if exist(fullfile(parentDir, currFileName), 'file')
           fileChecklistTable{iExp, iFile} = 1;
        end
    end
end

% Add expID column and save csv
fileChecklistTable = [table(expList, 'VariableNames', {'expID'}), fileChecklistTable];
writetable(fileChecklistTable, fullfile(saveDir, 'outputFileChecklist'));