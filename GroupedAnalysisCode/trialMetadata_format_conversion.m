parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';

expList = load_expList();
varCounts = [];
for iExp = 1:size(expList, 1)
    currFile = fullfile(parentDir, [expList.expID{iExp}, '_trialMetadata.csv']);
    disp(currFile);
    if exist(currFile, 'file')

        trialMd = readtable(currFile, 'delimiter', ',');
        if numel(trialMd.Properties.VariableNames) == 11
            shutoffVols = trialMd.pmtShutoffVols;
            trialMd.pmtShutoffVols = repmat({[]}, size(trialMd, 1), 1);
            for iTrial = 1:numel(shutoffVols)
                if ~isnan(shutoffVols(iTrial))
                    trialMd.pmtShutoffVols{iTrial} = {shutoffVols(iTrial)};
                end
            end
        elseif numel(trialMd.Properties.VariableNames) > 11
            baseData = trialMd(:, 1:10);
            shutoffArray = table2array(trialMd(:, 11:end));
            shutoffArray(isnan(shutoffArray)) = 0;
            shutoffVols = {};
            for iRow = 1:size(shutoffArray, 1)
                shutoffVols{iRow, 1} = find(shutoffArray(iRow, :));
            end
            trialMd = baseData;
            trialMd.pmtShutoffVols = shutoffVols;
        end
        trialMetadata = trialMd;
        save(fullfile(parentDir, [expList.expID{iExp}, '_trialMetadata.mat']), 'trialMetadata', ...
            '-v7');
        try
            movefile(currFile, fullfile(parentDir, 'trialMetadataBackup', [expList.expID{iExp}, ...
                '_trialMetadata.csv']));
        catch
            copyfile(currFile, fullfile(parentDir, 'trialMetadataBackup', [expList.expID{iExp}, ...
                '_trialMetadata.csv']));
        end
    end
end


%%
test2 = table2array(test(:, 11:end));
test2(isnan(test2)) = 0;

newData = {};
for iRow = 1:size(test2, 1)
    newData{iRow, 1} = find(test2(iRow, :));
end

test3 = test(:, 1:10);
test3.pmtShutoffVols = newData;
    
writetable(test3, 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\test.csv')
    
    
    
    
    