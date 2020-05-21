
% Just attempting to load and concatenate various tables for all experiments to make sure 
% they have identical fields

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';
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

% % Test each .csv file
% testTables = [];
% for iFile = 1:numel(uniqueFiles)
%     if contains(uniqueFiles{iFile}, '.csv')
%         disp(uniqueFiles{iFile})
%         testTables{end + 1} = [];
%         for iExp = 1:numel(expList)
%             currFile = fullfile(parentDir, [expList{iExp}, '_', uniqueFiles{iFile}]);
%             if exist(currFile, 'file')
%                 currData = readtable(currFile, 'delimiter', ',');
%                 testTables{end} = [testTables{end}; currData(1, :)];
%             end
%         end
%     end
% end

%%

uniqueFiles = uniqueFiles(contains(uniqueFiles, '.mat'));

blockMd = [];
ft = [];
fullRefImg = [];
roiData = [];
panelsMd = [];
subRoiMd = [];


for iExp = 1:numel(expList)
    
%     % blockMetadata.mat
%     currFile = fullfile(parentDir, [expList{iExp}, '_blockMetadata.mat']);
%     if exist(currFile, 'file')
%         load(currFile, 'blockMetadata');
%         blockMd = [blockMd; blockMetadata(1, :)];
%     end
    

%     % ficTracData.mat
%     currFile = fullfile(parentDir, [expList{iExp}, '_ficTracData.mat']);
%     if exist(currFile, 'file')
%         disp(currFile)
%         load(currFile, 'ftData');
%         fieldNames = fieldnames(ftData);
%         for iField = 1:numel(fieldNames)
%             if size(ftData(1).(fieldNames{iField}), 1) < size(ftData(1).(fieldNames{iField}), 2)
%                 for iRow = 1:numel(ftData)
%                     ftData(iRow).(fieldNames{iField}) = ftData(iRow).(fieldNames{iField})';
%                 end
%             end
%         end
%         ftData = struct2table(ftData, 'asarray', 1);
%         varNames = ftData.Properties.VariableNames;
%         if ~strcmp(varNames{1}, 'trialNum')
%             ftData = [table((1:size(ftData, 1))', 'VariableNames', {'trialNum'}), ftData];
%         end
%         save(currFile, 'ftData');
%     end
%     currFile = fullfile(parentDir, [expList{iExp}, '_ficTracData.mat']);
%     if exist(currFile, 'file')
%         disp(currFile)
%         load(currFile, 'ftData');
%         ft = [ft; ftData(1, :)];
%     end
    
%     % fullExpRefImages.mat
%     currFile = fullfile(parentDir, [expList{iExp}, '_fullExpRefImages.mat']);
%     if exist(currFile, 'file')
%         disp(currFile);
%         load(currFile, 'fullExpRefImages');
%     end
    
    % panelsMetadata.mat
    
%     refImages.mat
%         currFile = fullfile(parentDir, [expList{iExp}, '_refImages.mat']);
%         if exist(currFile, 'file')
%             disp(currFile);
%             load(currFile, 'refImages');
%             disp(size(refImages));
%         end
    % roiData.mat
    
    % subRoiMetadata.mat

end





