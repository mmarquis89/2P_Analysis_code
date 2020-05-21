
% Load exp list
expList = load_expList('groupName', 'gaplessAcq');
% expList = expList(1:17,:)
% expList = expList(~cellfun(@isempty, regexp(expList.expID, '2018.*', 'match', 'once')), 1);
% expList = load_expList();
% expList = expList((~cellfun(@isempty, (regexp(expList.expID, '201[78].*')))),:)

summaryStats = [];
for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID);
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    parentDir = find_parent_dir(currExpID);
    newRow = table({currExpID}, 'VariableNames', {'expID'});
    if exist(fullfile(parentDir, expDirName, 'volCounts.mat'), 'file')
        load(fullfile(parentDir, expDirName, 'volCounts.mat'), 'volCounts');
        load(fullfile(parentDir, expDirName, 'summaryStats.mat'), 'medVals', 'minVals', ...
                'minFrameVals')
        
%         newRow.minVals = {cell2mat(minVals')};
        newRow.medVals = {cell2mat(medVals')};
        newRow.volCounts = {volCounts};
        newRow.trialBounds = {cumsum(volCounts)};
        trialMedVals = [];
        trialStarts = [1, newRow.trialBounds{:}(1:end-1)+1];
        for iTrial = 1:numel(newRow.trialBounds{:})
            trialMedVals(iTrial) = median(newRow.medVals{:}( ...
                    trialStarts(iTrial):newRow.trialBounds{:}(iTrial)));
        end
        newRow.trialMedVals = {trialMedVals};
        
    else
        newRow.medVals = {[]};
        newRow.volCounts = {[]};
        newRow.trialBounds = {[]};
        newRow.trialMedVals = {[]};
    end
    summaryStats = [summaryStats; newRow];
end
%%

test = summaryStats.trialMedVals;
test2 = cellfun(@(x) max([0, (max(x) - min(x))]), test);

test3 = table(summaryStats.expID, test2, 'VariableNames', {'expID', 'maxOffset'});

%%
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';

expList = load_expList();

nRows = [];
nCols = [];
for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID)
    roiFile = fullfile(parentDir, [currExpID, '_roiData.mat']);
    if exist(roiFile, 'file')
        
        % Load ROI data for current experiment
        load(roiFile, 'roiData');

        nRows(iExp) = size(roiData.rawFl{1}, 1);
        nCols(iExp) = size(roiData.rawFl{1}, 2);
    else
        nRows(iExp) = nan;
        nCols(iExp) = nan;
    end%if
end

newTable = [expList, table(nRows', nCols', 'VariableNames', {'nRows', 'nCols'})];












