parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';

expList = load_expList();

% Load FicTrac and ROI data for all experiments
allExpData = [];
for iExp = 1:size(expList, 1)
    
    currExpID = expList.expID{iExp};
    disp(currExpID);
    
    roiFileName = [currExpID, '_roiData.mat'];
    roiFile = fullfile(parentDir, roiFileName);
    ftFileName = [currExpID, '_ficTracData.mat'];
    ftFile = fullfile(parentDir, ftFileName);
    
    expMdFile = fullfile(parentDir, [currExpID, '_expMetadata.csv']);
    trialMdFile = fullfile(parentDir, [currExpID, '_trialMetadata.mat']);

    
    if exist(roiFile, 'file') && exist(ftFile, 'file') && exist(expMdFile, 'file')
       
        
        load(roiFile, 'roiData');
        load(ftFile, 'ftData');
        expMd = readtable(expMdFile, 'delimiter', ',');
        load(trialMdFile, 'trialMetadata');
        trialMd = trialMetadata;
        %         % Fix FicTrac data if I accidentally saved it as a struct
%         if isstruct(ftData)
%             copyfile(ftFile, fullfile(parentDir, 'ftBackup', ftFileName));
%             
%             ftData = struct2table(ftData, 'asarray', 1);
%             if strcmp(ftData.Properties.VariableNames{1}, 'trialNum')
%                ftData = ftData(:, 2:end); 
%             end
%             for iRow = 1:size(ftData, 1)
%                 sz = size(ftData.frameTimes{iRow});
%                 if sz(2) > sz(1)
%                    ftData.frameTimes{iRow} = ftData.frameTimes{iRow}'; 
%                 end
%             end
%             expID = repmat({currExpID}, size(ftData, 1), 1);
%             trialNum = (1:size(ftData, 1))';
%             ftData = [table(expID, trialNum, 'VariableNames', {'expID', 'trialNum'}), ftData];
%             save(ftFile, 'ftData');
%         else
%             changed = 1;
%             for iRow = 1:size(ftData, 1)
%                 sz = size(ftData.frameTimes{iRow});
%                 if sz(2) > sz(1)
%                     ftData.frameTimes{iRow} = ftData.frameTimes{iRow}';
%                     changed = 1;
%                 end
%             end
%             if changed
%                 copyfile(ftFile, fullfile(parentDir, 'ftBackup', ftFileName));
%                 save(ftFile, 'ftData');
%             end
%         end
%         
%         % Fix ROI data if I accidentally included dF/F fields
%         if any(strcmp(roiData.Properties.VariableNames, 'dffData'))
%             
%             copyfile(roiFile, fullfile(parentDir, 'roiDataBackup', roiFileName));
%             
%             fixedRoiData = roiData(:, {'expID', 'trialNum', 'roiName', 'subROIs', 'rawFl'});
%             
%             % Replace trial-based dF/F field with just the trial-based baseline values
%             trialBaseline = [];
%             for iRow = 1:size(roiData, 1)
%                 roiDataSorted = sort(roiData.rawFl{iRow});
%                 trialBaseline(iRow, 1) = median(roiDataSorted(1:round(numel(roiDataSorted) * 0.05)));
%             end
%             fixedRoiData.trialBaseline = trialBaseline;
%             
%             % Replace expDff field with the full-exp baseline
%             roiList = unique(fixedRoiData.roiName);
%             baselineVals = [];
%             for iROI = 1:size(roiList, 1)
%                 currRoiData = fixedRoiData(strcmp(roiData.roiName, roiList{iROI}), :);
%                 roiDataSort = sort(cell2mat(currRoiData.rawFl));
%                 baselineVals(iROI) = median(roiDataSort(1:round(numel(roiDataSort) * 0.05)));
%             end
%             baselineTable = [table(roiList, 'variablenames', {'roiName'}), table(baselineVals', ...
%                     'variableNames', {'expBaseline'})];
%             roiData = innerjoin(fixedRoiData, baselineTable);            
%             
%             save(roiFile, 'roiData');
%         end
        
        allExpData = [allExpData; inner_join(expMd, trialMd, roiData, ftData)];
    end    
end

%%

dt = DataTable(allExpData);
filterDefs = struct();
filterDefs.roiName = 'TypeF';
filterDefs.expID = '20190226-3'%'20181020-1';%

smWin = 7;

dt = dt.initialize_filters(filterDefs);

%

currData = dt.subset(:, {'expID', 'trialNum', 'roiName', 'volumeRate', 'nVolumes', 'trialDuration', ...
        'originalTrialCount', 'rawFl', 'trialBaseline', 'expBaseline', 'moveSpeed', 'yawSpeed', ...
        'fwSpeed', 'frameTimes', 'pmtShutoffVols'});

% Calculate volume times
volTimes = calc_volTimes(currData.nVolumes, currData.volumeRate, currData.trialDuration, ...
        currData.originalTrialCount);
    
% Get smoothed move and yaw speed data
smMoveSpeed = smoothdata(currData.moveSpeed{1}, 1, 'gaussian', smWin);

% Convert move speed from rad/frame to mm/sec 
smMoveSpeed = smMoveSpeed * 25 * 4.5;

% Downsample smoothed data to match volTimes  
speedDataVols = [];
for iVol = 1:numel(volTimes)
    speedDataVols(iVol) = smMoveSpeed(argmin(abs(currData.frameTimes{:} - volTimes(iVol))));
end
speedDataVols = smoothdata(speedDataVols', 1, 'gaussian', smWin);

% Get smoothed fl data
smFl = smoothdata(currData.rawFl{1}, 1, 'gaussian', smWin);

% Drop any nan volumes
speedDataVols(currData.pmtShutoffVols{1}) = nan;
smFl(currData.pmtShutoffVols{1}) = nan;

% Apply 3-volume lag
speedDataVols = speedDataVols(1:end - 3);
smFl = smFl(4:end);

% Overlay normalized data
figure(1);clf;hold on
set(gcf, 'Color', [1 1 1]);
plot(speedDataVols ./ max(speedDataVols));
plot(smFl ./ max(smFl));

% Scatterplot
figure(2);clf;
set(gcf, 'Color', [1 1 1]);
plot(smFl, speedDataVols, '.')
xlabel('rawFl'); ylabel('moveSpeed')

% 2D histogram
figure(3);clf;
set(gcf, 'Color', [1 1 1]);
speedDataVols(smFl < 100) = nan;
smFl(smFl < 100) = nan;
histogram2(smFl, speedDataVols, 45, 'displaystyle', 'tile')









%%


corrcoef




