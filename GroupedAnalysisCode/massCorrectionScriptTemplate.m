parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';

expList = load_expList();

if ~isdir(fullfile(parentDir, 'roiDataBackup'))
    mkdir(fullfile(parentDir, 'roiDataBackup'))
end

for iExp = 1:size(expList, 1)
    
    currExpID = expList.expID{iExp};
    roiDataFileName = fullfile(parentDir, [currExpID, '_roiData.mat']);
    if exist(roiDataFileName, 'file')
        copyfile(roiDataFileName, fullfile(parentDir, 'roiDataBackup', ...
                [currExpID, '_roiData.mat']));
        load(roiDataFileName, 'roiData');
        
        % Remove dffData and expDffData fields
        roiData = roiData(:, 1:end-2);
        
        % Save trial baseline values for future dF/F calculation
        trialBaseline = [];
        for iRow = 1:size(roiData, 1)
            roiDataSorted = sort(roiData.rawFl{iRow});
            trialBaseline(iRow) = median(roiDataSorted(1:round(numel(roiDataSorted) * 0.05)));
        end
        roiData.trialBaseline = trialBaseline';
        
        % Save full experiment baseline values
        expBaseline = nan(size(roiData, 1), 1);
        roiList = unique(roiData.roiName);
        for iRoi = 1:numel(roiList)
            currRoiData = roiData(strcmp(roiData.roiName, roiList{iRoi}), :);
            flDataSort = sort(cell2mat(currRoiData.rawFl));
            currRoiBaseline = repmat(median(flDataSort(1:round(numel(flDataSort) * 0.05))), ...
                    size(currRoiData, 1), 1);
            expBaseline(strcmp(roiData.roiName, roiList{iRoi})) = currRoiBaseline;
        end       
        roiData.expBaseline = expBaseline;
        
        
        save(roiDataFileName, 'roiData', '-v7');
    end    
    
end

%% 

dffCalc = @(x, y) (x - y)./y;
test = roiData(1,:)

tempDff = dffCalc(test.rawFl{:}, test.trialBaseline);
tempExpDff = dffCalc(test.rawFl{:}, test.expBaseline);

figure(1);clf;hold on
plot(tempDff);
plot(test.dffData{:});

figure(2);clf;hold on;
plot(tempExpDff);
plot(test.expDffData{:});



