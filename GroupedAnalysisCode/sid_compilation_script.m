parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2019 Jan-Feb';

expDate = '2019_02_16_exp_1';
sidNums = [0 1];

expDir = fullfile(parentDir, expDate);
saveDir = fullfile(expDir, 'sid_master');
if ~isdir(saveDir)
    mkdir(saveDir)
end

% Load files for each sid
aD = {}; volAvgData = {}; ROIdata = {}; ROImd = {};
shortTrialCounts = []; blockCounts = [];
for iSid = 1:numel(sidNums)
    currDir = fullfile(expDir, ['sid_', num2str(sidNums(iSid))]);
    load(fullfile(currDir, 'analysisMetadata.mat'), 'analysisMetadata');
    load(fullfile(currDir, ['sid_', num2str(sidNums(iSid)), '_volAvgSessionData.mat']), ...
            'volAvgSessionData');
    load(fullfile(currDir, 'ROI_Data_Avg.mat'), 'ROIDataAvg');
    load(fullfile(currDir, 'ROI_metadata.mat'), 'ROImetadata');
    disp([analysisMetadata.volumeRate analysisMetadata.nPlanes ...
            analysisMetadata.nVolumes analysisMetadata.trialDuration  numel(ROImetadata) ...
            numel(analysisMetadata.blockData) analysisMetadata.nTrials]);
        
    aD{iSid} = analysisMetadata;
    shortTrialCounts(iSid) = analysisMetadata.nTrials;
    blockCounts(iSid) = numel(analysisMetadata.blockData);
    volAvgData{iSid} = volAvgSessionData;
    ROIdata{iSid} = ROIDataAvg;
    ROImd{iSid} = ROImetadata; 
end

% Create a table to track the number of short trials that belong to each session
sidTrialCounts = table(sidNums', shortTrialCounts', blockCounts', ...
        'VariableNames', {'sidNum', 'nShortTrials', 'nBlocks'});

    
% ---- Now the concatenation steps -----

% Copy over the required info that is shared between sids
aD_master = struct();
aD_master.expDate = aD{1}.expDate;
aD_master.volumeRate = aD{1}.volumeRate;
aD_master.nPlanes = aD{1}.nPlanes;
aD_master.trialDuration = aD{1}.trialDuration;
aD_master.nVolumes = aD{1}.nVolumes;

for iSid = 1:numel(sidNums)
    
    % Just copy everything over for the first sid
    if iSid == 1
        aD_master.blockData = aD{1}.blockData;
        aD_master.flowArr = aD{1}.flowArr;
        aD_master.ftData = aD{1}.ftData;
        aD_master.trialAnnotations = aD{1}.trialAnnotations;
        aD_master.nTrials = aD{1}.nTrials;
        volAvgSessionData_master = volAvgData{1};
        ROIDataAvg_master = ROIdata{1};
        ROImetadata_master = ROImd{1};
    else
    
        % Concatenate the rest onto the first one
        aD_master.nTrials = aD_master.nTrials + aD{iSid}.nTrials;
        
        % blockData
        aD_master.blockData = [aD_master.blockData, aD{iSid}.blockData];
        
        % ftData/flow data
        aD_master.flowArr = [aD_master.flowArr, aD{iSid}.flowArr];
        ftFieldNames = fieldnames(aD_master.ftData);
        for iField = 1:numel(ftFieldNames)
            currField = ftFieldNames{iField};
            sz = size(aD_master.ftData.(currField));
            if strcmp(currField, 'badFtTrials')
                aD_master.ftData.badFtTrials = ...
                    [aD_master.ftData.badFtTrials, aD_master.ftData.badFtTrials + aD{iSid}.nTrials];
            elseif sz(1) == 1
                aD_master.ftData.(currField) = [aD_master.ftData.(currField), ...
                        aD{iSid}.ftData.(currField)];
            else
                aD_master.ftData.(currField) = ...
                    cat(numel(sz), aD_master.ftData.(currField), aD{iSid}.ftData.(currField));
            end
        end
        
        % annotation data
        aD_master.trialAnnotations = [aD_master.trialAnnotations; aD{iSid}.trialAnnotations];
        
        % volAvgSessionData
        volAvgSessionData_master = cat(4, volAvgSessionData_master, volAvgData{iSid});
        
        % ROIDataAvg
        ROIDataAvg_master = cat(2, ROIDataAvg_master, ROIdata{iSid});
        
    end%if
end%iSid

% ----- Rename and save in master directory -----
analysisMetadata = aD_master;
volAvgSessionData = volAvgSessionData_master;
ROIDataAvg = ROIDataAvg_master;
ROImetadata = ROImd;
save(fullfile(saveDir, 'analysisMetadata.mat'), 'analysisMetadata');
save(fullfile(saveDir, 'volAvgSessionData.mat'), 'volAvgSessionData');
save(fullfile(saveDir, 'ROI_Data_Avg.mat'), 'ROIDataAvg');
for iSid = 1:numel(ROImd)
    ROImetadata = ROImd{iSid};
    fileName = ['sid_', num2str(sidNums(iSid)), '_ROI_metadata.mat'];
    save(fullfile(saveDir, fileName), 'ROImetadata');
end
writetable(sidTrialCounts, fullfile(saveDir, 'sidTrialCounts.csv'));





