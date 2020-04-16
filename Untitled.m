

%% Reformat expIDs and names from my manual Excel spreadsheet
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
expList = readtable(fullfile(parentDir, 'oldExpList_gaplessAcq.csv'), 'delimiter', ',');

years = regexp(expList.year, '..(?= \()', 'match', 'once');
expNums = regexp(expList.year, '.(?=\))', 'match', 'once');


newExpList = [];
for iExp = 1:numel(years)
    expID = ['20', years{iExp}, pad(num2str(expList.month(iExp)), 2, 'left', '0'), ...
            pad(num2str(expList.day(iExp)), 2, 'left', '0'), '-', expNums{iExp}];
    newRow = table({expID}, expList.expName(iExp), 'VariableNames', {'expID', 'expName'});
    newExpList = [newExpList; newRow];        
end

writetable(newExpList, fullfile(parentDir, 'oldExpList_gaplessAcq.csv'))

%% Loop through each experiment and double check the format of some key data

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
expList = readtable(fullfile(parentDir, 'oldExpList_gaplessAcq.csv'), 'delimiter', ',');

% expList = expList(20:end, :); % temporarily skip older exps that aren't available on local HD
% expList = expList([20, 29:size(expList, 1)], :); % also skipping multi-sid (nerve transection) expts
%  expList = expList(21:28, :); % focus specifically on nerve transection expts

dataParentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';
oldExpParentDirs = dir(fullfile(dataParentDir, '2019 *'));
expInfo = struct();

for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    expInfo(iExp).expID = currExpID;
    disp(currExpID);
    
    % Convert expID into name of experiment dir
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    
    % Find experiment directory
    if ~isempty(dir(fullfile(dataParentDir, expDirName)))
        expDir = fullfile(dataParentDir, expDirName);
    elseif ~isempty(dir(fullfile(dataParentDir, '2018', expDirName)))
        expDir = fullfile(dataParentDir, '2018', expDirName);
    else
        for iDir = 1:numel(oldExpParentDirs)
            currDir = fullfile(oldExpParentDirs(iDir).folder, oldExpParentDirs(iDir).name);
            expDir = fullfile(currDir, expDirName);
            if ~isempty(dir(expDir))
                break
            end
        end        
    end
    
    % Identify number of session folders
    sidDirs = dir(fullfile(expDir, 'sid*'));
    sidDirs = sidDirs([sidDirs.isdir]);
    expInfo(iExp).nSids = numel(sidDirs);
    sidDirNames = {sidDirs.name};
    expInfo(iExp).sidNums = cellfun(@(x) x(end), sidDirNames);
    
    % Load analysis data struct (from first session only if there's more than one)
    analysisDir = fullfile(sidDirs(1).folder, sidDirs(1).name);
    load(fullfile(analysisDir, 'analysisMetadata.mat'), 'analysisMetadata');
    aD = analysisMetadata;
    
    % Compare number of trials in block data vs. full exp data
    expInfo(iExp).nTrials = aD.nTrials;
    expInfo(iExp).nTrialsDaq = sum([aD.blockData.nTrials]);
    
    % Record number of output channels in DAQ data 
    expInfo(iExp).nDaqOutputs = size(aD.blockData(1).outputData, 2);
    
    % Record number of input and output daq samples for first block data
    expInfo(iExp).nDaqSamplesOutput = size(aD.blockData(1).outputData, 1);
    expInfo(iExp).nDaqSamplesInput = size(aD.blockData(1).inputData, 1);
    
    % Record list of stim types 
    expInfo(iExp).stimTypes = aD.stimTypes;
    
    % Record total number of fields in analysis metadata structure
    expInfo(iExp).nStructFields = numel(fieldnames(aD));
    
end 











