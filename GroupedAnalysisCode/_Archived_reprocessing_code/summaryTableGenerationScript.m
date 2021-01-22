
saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
parentDir = fullfile(saveDir, 'all_experiments');

fullExpList = readtable(fullfile(saveDir, 'fullExpList.csv'), 'delimiter', ',');
oneNoteMd = readtable(fullfile(saveDir, 'OneNoteMetadata.txt'));

expMdFiles = dir(fullfile(parentDir, '*expMetadata.csv'));
allExpMd = [];
for iFile = 1:numel(expMdFiles)
   currMd = readtable(fullfile(parentDir, expMdFiles(iFile).name), 'delimiter', ',');
   allExpMd = [allExpMd; currMd];    
end

fullExpMdTable = innerjoin(innerjoin(fullExpList, allExpMd), oneNoteMd);

%

expList = fullExpMdTable.expID;
optoStim = [];
odorStim = [];
panelsFlashStim = [];
soundStim = [];
laserStim = [];
for iExp = 1:numel(expList)
    
    currExpID = expList{iExp};
    eventDataFiles = dir(fullfile(parentDir, [currExpID, '_event_data_*.csv']));
    if any(~cellfun(@isempty, regexp({eventDataFiles.name}, 'optostim')))
        optoStim(iExp) = 1;
    else 
        optoStim(iExp) = 0;
    end
    if any(~cellfun(@isempty, regexp({eventDataFiles.name}, 'odor')))
        odorStim(iExp) = 1;
    else 
        odorStim(iExp) = 0;
    end
    if any(~cellfun(@isempty, regexp({eventDataFiles.name}, 'panelsflash')))
        panelsFlashStim(iExp) = 1;
    else 
        panelsFlashStim(iExp) = 0;
    end
    if any(~cellfun(@isempty, regexp({eventDataFiles.name}, 'soundstim')))
        soundStim(iExp) = 1;
    else 
        soundStim(iExp) = 0;
    end
    if any(~cellfun(@isempty, regexp({eventDataFiles.name}, 'ir-laser')))
        laserStim(iExp) = 1;
    else
        laserStim(iExp) = 0;
    end
end

fullExpMdTable.optoStim = optoStim';
fullExpMdTable.odorStim = odorStim';
fullExpMdTable.panelsFlashStim = panelsFlashStim';
fullExpMdTable.soundStim = soundStim';
fullExpMdTable.IR_laserStim = laserStim';

% 

odorNames = {};
for iExp = 1:numel(expList)
    currExpID = expList{iExp};
    if fullExpMdTable.odorStim(iExp)
        odorEventData = readtable(fullfile(parentDir, [currExpID, '_event_data_odor.csv']), ...
                'delimiter', ',');
        uniqueNames = unique(odorEventData.odorName);
        odorNames{iExp, 1} = uniqueNames{1};
        if numel(uniqueNames) > 1
           for iName = 2:numel(uniqueNames)
               odorNames{iExp, 1} = [odorNames{iExp, 1}, ', ', uniqueNames{iName}];
           end
        end
    else
        odorNames{iExp, 1} = [];
    end
end

fullExpMdTable.odorNames = odorNames;

%

roiNames = {};
for iExp = 1:numel(expList)
    currExpID = expList{iExp};
    if exist(fullfile(parentDir, [currExpID, '_roiData.mat']), 'file')
        load(fullfile(parentDir, [currExpID, '_roiData.mat']), 'roiData');
        uniqueNames = unique(roiData.roiName);
        roiNames{iExp, 1} = uniqueNames{1};
        if numel(uniqueNames) > 1
           for iName = 2:numel(uniqueNames)
               roiNames{iExp, 1} = [roiNames{iExp, 1}, ', ', uniqueNames{iName}];
           end
        end
    else
        roiNames{iExp, 1} = [];
    end
end

fullExpMdTable.roiNames = roiNames;

%%

for i=1:numel(roiNames)
   if isempty(roiNames{i})
       roiNames{i} = '';
   end
end

roiNameTable = table(expList, 'VariableNames', {'expID'});

roiNameTable.Background = ~cellfun(@isempty, regexp(roiNames, 'Background'));
roiNameTable.ANT_L = ~cellfun(@isempty, regexp(roiNames, 'ANT-L'));
roiNameTable.ANT_R = ~cellfun(@isempty, regexp(roiNames, 'ANT-R'));
roiNameTable.TypeD_L = ~cellfun(@isempty, regexp(roiNames, 'TypeD-L'));
roiNameTable.TypeD_R = ~cellfun(@isempty, regexp(roiNames, 'TypeD-R'));
roiNameTable.TypeB_L = ~cellfun(@isempty, regexp(roiNames, 'TypeB-L|LH-L'));
roiNameTable.TypeB_R = ~cellfun(@isempty, regexp(roiNames, 'TypeB-R|LH-R'));
roiNameTable.TypeC_L = ~cellfun(@isempty, regexp(roiNames, 'TypeC-L'));
roiNameTable.TypeC_R = ~cellfun(@isempty, regexp(roiNames, 'TypeC-R'));
roiNameTable.TypeF_L = ~cellfun(@isempty, regexp(roiNames, 'TypeF-L'));
roiNameTable.TypeF_R = ~cellfun(@isempty, regexp(roiNames, 'TypeF-R'));
roiNameTable.VLP_AMMC_L = ~cellfun(@isempty, regexp(roiNames, 'VLP-AMMC-L'));
roiNameTable.VLP_AMMC_R = ~cellfun(@isempty, regexp(roiNames, 'VLP-AMMC-R'));
roiNameTable.SMP_DN =~cellfun(@isempty, regexp(roiNames, 'SMP-DN'));
roiNameTable.PB = ~cellfun(@isempty, regexp(roiNames, 'PB'));

roiNameTable.uni_ANT = ~cellfun(@isempty, regexp(roiNames, 'ANT'));
roiNameTable.bi_ANT = roiNameTable.ANT_L & roiNameTable.ANT_R;

roiNameTable.uni_TypeD = ~cellfun(@isempty, regexp(roiNames, 'TypeD'));
roiNameTable.bi_TypeD = roiNameTable.TypeD_L & roiNameTable.TypeD_R;

roiNameTable.uni_TypeB = ~cellfun(@isempty, regexp(roiNames, 'TypeB|LH'));
roiNameTable.bi_TypeB = roiNameTable.TypeB_L & roiNameTable.TypeB_R;

roiNameTable.uni_TypeC = ~cellfun(@isempty, regexp(roiNames, 'TypeC'));
roiNameTable.bi_TypeC = roiNameTable.TypeC_L & roiNameTable.TypeC_R;

roiNameTable.uni_TypeF = ~cellfun(@isempty, regexp(roiNames, 'TypeF'));
roiNameTable.bi_TypeF = roiNameTable.TypeF_L & roiNameTable.TypeF_R;

roiNameTable.uni_VLP_AMMC = roiNameTable.VLP_AMMC_L | roiNameTable.VLP_AMMC_R;
roiNameTable.bi_VLP_AMMC = roiNameTable.VLP_AMMC_L & roiNameTable.VLP_AMMC_R;

writetable(roiNameTable, fullfile(saveDir, 'roiNameTable.csv'))










