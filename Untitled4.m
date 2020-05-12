
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

test3 = table(summaryStats.expID, test2, 'VariableNames', {'expID', 'maxOffset'})


%%

for iExp = 1:size(summaryStats, 1)
    f = figure(iExp); clf; hold on
    f.Color = [1 1 1];
    f.Position = [50, 250, 1500, 700];
    %    plot(smoothdata(summaryStats.minVals{iExp}, 'gaussian', 5));
    %    plot(smoothdata(summaryStats.medVals{iExp}, 'gaussian', 5));
    plot(summaryStats.trialMedVals{iExp})
    yL = ylim();
%     trialBounds = cumsum(summaryStats.volCounts{iExp});
%     for iTrial = 1:numel(trialBounds)
%         plot([trialBounds(iTrial), trialBounds(iTrial)], yL, '--', 'color', 'm')
%     end
%     ax = gca;
%     ax.XTick = trialBounds;
%     ax.XTickLabel = 1:numel(trialBounds);
    title(summaryStats.expID{iExp});
end



%%

% Load exp list
expList = load_expList();
% expList = expList(~cellfun(@isempty, regexp(expList.expID, '2018.*', 'match', 'once')), 1);

summaryStats = [];
for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID);
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    parentDir = find_parent_dir(currExpID);
    
    if exist(fullfile(parentDir, expDirName, 'volCounts.mat'), 'file')
        load(fullfile(parentDir, expDirName, 'volCounts.mat'), 'volCounts');
        load(fullfile(parentDir, expDirName, 'summaryStats.mat'), 'medVals', 'minVals', ...
                'minFrameVals')
        newRow = table({currExpID}, 'VariableNames', {'expID'});
        newRow.minVals = {cell2mat(minVals')};
        newRow.medVals = {cell2mat(medVals')};
        newRow.volCounts = {volCounts};
        summaryStats = [summaryStats; newRow];
    end
end

%%

for iExp = 1:size(summaryStats, 1)    
   f = figure(iExp); clf; hold on
   f.Color = [1 1 1];
   f.Position = [50, 250, 1500, 700];
   plot(smoothdata(summaryStats.minVals{iExp}, 'gaussian', 5));
   plot(smoothdata(summaryStats.medVals{iExp}, 'gaussian', 5));
   yL = ylim();
   trialBounds = cumsum(summaryStats.volCounts{iExp});
   for iTrial = 1:numel(trialBounds)
       plot([trialBounds(iTrial), trialBounds(iTrial)], yL, '--', 'color', 'm')
   end
   ax = gca;
   ax.XTick = trialBounds
   ax.XTickLabel = 1:numel(trialBounds)
   title(summaryStats.expID{iExp});
end

%%
expList = load_expList('groupName', 'singleTrialAcq_preROIs');
% expList = expList(4:end, :);
for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID);
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
    parentDir = find_parent_dir(currExpID);
    refImgFiles = dir(fullfile(parentDir, expDirName, 'refImages_reg_block*.mat'));
    if ~isempty(refImgFiles)
        for iFile = 1:numel(refImgFiles)
            load(fullfile(parentDir, expDirName, refImgFiles(iFile).name), 'refImages');
            refImages = squeeze(mean(refImages, 5));
            save(fullfile(parentDir, expDirName, refImgFiles(iFile).name), 'refImages');
        end
    end
end


















