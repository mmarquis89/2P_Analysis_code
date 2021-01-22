saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\'; 
% expList = load_expList('groupname', 'gaplessAcq');
% 
% % Skipping multi-sid (nerve transection) expts because I already did them
% expList = expList([5 11], :); 

expList = table({'20190315-3'}, 'variablenames', {'expID'}); 
expList.expName{1} = 'D-ANT_6s';

expNum =1;

% ----- AUTOMATIC PROCESSING AND CONVERSION -----
try

% Convert expID to old-format expDir name
currExpID = expList.expID{expNum};
disp(currExpID);
expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];

% Find experiment directory
expParentDir = find_parent_dir(currExpID);
sidDir = dir(fullfile(expParentDir, expDirName, 'sid_*'));
sidDir = sidDir([sidDir.isdir]);
if numel(sidDir) == 1
    expDir = fullfile(expParentDir, expDirName, sidDir.name); % If there's only one sid directory
elseif numel(sidDir) > 1
    expDir = fullfile(expParentDir, expDirName, 'sid_master'); % If I've already manually combined sids
else
    error('Could not find experiment directory');
end
expName = expList.expName{expNum};

% Load analysisMetadata file
load(fullfile(expDir, 'analysisMetadata.mat'), 'analysisMetadata');
aD = analysisMetadata;

FRAME_RATE = 25; % Behavior video (and therefore FicTrac) frame rate was a constant

% ----- CREATE EXP METADATA TABLE -----

expID = regexprep(aD.expDate, {'_', 'exp'}, {'', '-'});
expMd = table({expID}, 'VariableNames', {'expID'});
expMd.expName = {expName};
expMd.daqSampRate = 4000;        % this is after 10x downsampling at the time of the experiment, and will be changed later in the script to downsample even more if saving DAQ output data
expMd.panelsDisplayRate = 50;    % Didn't have panels but just keeping the number consistent
expMd.volumeRate = aD.volumeRate;
expMd.nPlanes = aD.nPlanes; 
expMd.nTrials = numel(aD.blockData); % Redefining "trial" as an individual continuous acquisition

% ----- CREATE TRIAL METADATA TABLE -----
trialMd = [];
for iTrial = 1:expMd.nTrials
    newRow = table({expID}, 'VariableNames', {'expID'});
    newRow.trialNum = iTrial;
    newRow.trialDuration = aD.trialDuration * aD.blockData(iTrial).nTrials;
    newRow.nVolumes = aD.nVolumes * aD.blockData(iTrial).nTrials;
    newRow.nDaqSamples = size(aD.blockData(iTrial).outputData, 1);
    newRow.nPanelsFrames = expMd.panelsDisplayRate * newRow.trialDuration;
    newRow.usingOptoStim = 0;
    newRow.usingPanels = 0;
    newRow.using2P = 1;
    newRow.originalTrialCount = aD.blockData(iTrial).nTrials;
    newRow.pmtShutoffVols = {[]};
    trialMd = [trialMd; newRow];
end

% ----- DAQ OUTPUT DATA (DOWNSAMPLED FURTHER TO 100 Hz) -----
trialNum = [];
optoStimCommand = [];
speakerCommand = [];
odorACommand = [];
odorBCommand = [];
for iTrial = 1:expMd.nTrials
    currTrialOutputData = aD.blockData(iTrial).outputData;
    
    % 'Next file trigger' and 'acquisition start trigger' shared the same channel until 1/24/19,
    % when I split them from the first column into the first and second columns, so if the exp 
    % predates that insert a column of 'nan' on the left side of the array.
    if size(currTrialOutputData, 2) == 7
        currTrialOutputData = [nan(size(currTrialOutputData, 1), 1), currTrialOutputData];
    end
    
    % Opto stim column (3) was added starting 04/01/19...if exp predates that, insert a column 
    % of zeros where it should be    
    if size(currTrialOutputData, 2) == 8
        currTrialOutputData = [currTrialOutputData(:, 1:2), zeros(size(currTrialOutputData, 1), ...
                1), currTrialOutputData(:, 3:8)]; 
    end
    
    % After both corrections above column labels should be:
    % [acqStart, nextFile, optoStim, speaker, olfA, olfB, olfNO, trialAlignLED, cameraTrigger] 
    
    trialNum = [trialNum; ones(size(currTrialOutputData, 1), 1) * trialMd.trialNum(iTrial)];
    optoStimCommand = [optoStimCommand; currTrialOutputData(:, 3)];
    speakerCommand = [speakerCommand; currTrialOutputData(:, 4)];
    odorACommand = [odorACommand; currTrialOutputData(:, 5)];
    odorBCommand = [odorBCommand; currTrialOutputData(:, 6)];
        
end

% Smooth a copy of the opto stim command before downsampling so that it will never miss an opto
% stim pulse even if it's at a 1% duty cycle
dsFactor = expMd.daqSampRate / 100;
optoStimCommandSmoothed = smoothdata(optoStimCommand, 1, 'gaussian', 100);
dsOptoStimSmoothed = optoStimCommandSmoothed(1:dsFactor:end);

% Downsample
dsInds = 1:dsFactor:size(trialNum, 1);
daqOutputData = table(trialNum(dsInds), optoStimCommand(dsInds), speakerCommand(dsInds), ...
        odorACommand(dsInds), odorBCommand(dsInds), optoStimCommandSmoothed(dsInds), ...
        'VariableNames', {'trialNum', 'optoStim', 'speaker', 'odorA', 'odorB', 'optoStimSmoothed'});
expMd.daqSampRate = 100;


% ----- IDENTIFY BLOCK (AKA TRIAL) BOUNDARIES -----
nBlockTrials = [aD.blockData.nTrials];
blockStartTrials = 1 + [0, cumsum(nBlockTrials(1:end-1))];
blockEndTrials = cumsum(nBlockTrials);  


% ----- CREATE FICTRAC DATA STRUCTURE -----
ftData = struct();
for iBlock = 1:numel(nBlockTrials)
    currBlockTrials = blockStartTrials(iBlock):blockEndTrials(iBlock);
    ftData(iBlock).trialNum = iBlock;
    ftData(iBlock).intX = as_vector(squeeze(aD.ftData.intXY(:, 1, currBlockTrials)));
    ftData(iBlock).intY = as_vector(squeeze(aD.ftData.intXY(:, 2, currBlockTrials)));
    ftData(iBlock).intHD = as_vector(aD.ftData.intHD(:, currBlockTrials));
    ftData(iBlock).moveSpeed = as_vector(aD.ftData.moveSpeed(:, currBlockTrials));
    ftData(iBlock).intFwMove = as_vector(aD.ftData.intForwardMove(:, currBlockTrials));
    ftData(iBlock).intSideMove = as_vector(aD.ftData.intSideMove(:, currBlockTrials));
    ftData(iBlock).yawSpeed = as_vector(aD.ftData.yawSpeed(:, currBlockTrials));
    ftData(iBlock).fwSpeed = as_vector(aD.ftData.fwSpeed(:, currBlockTrials));
    ftData(iBlock).sideSpeed = as_vector(aD.ftData.sideSpeed(:, currBlockTrials));
    currBlockFrameCount = numel(ftData(iBlock).intX);
    ftData(iBlock).frameTimes = (1:currBlockFrameCount) / FRAME_RATE;
    ftData(iBlock).badVidFrames = isnan(ftData(iBlock).intX);
    ftData(iBlock).meanFlow = as_vector(aD.flowArr(:, currBlockTrials));
end


% ----- BEHAVIOR ANNOTATION EVENT DATA ----
quiescenceEvents = behaviorEvent('Quiescence');
isoMoveEvents = behaviorEvent('IsolatedMovement');
groomEvents = behaviorEvent('Grooming');
locEvents = behaviorEvent('Locomotion');
for iBlock = 1:numel(nBlockTrials)
    currBlockTrials = blockStartTrials(iBlock):blockEndTrials(iBlock);
    currAnnotData = aD.trialAnnotations(currBlockTrials, :)';
    currAnnotData(1, :) = currAnnotData(2, :); % currently the first frame is (wrongly) always quiescence
    currAnnotData = currAnnotData(:);
    frameTimes = (1:numel(currAnnotData)) / FRAME_RATE;
    
    % Separate according to behavior type
    qFrames = currAnnotData == 0;
    iFrames = currAnnotData == 1;
    gFrames = currAnnotData == 2;
    lFrames = currAnnotData == 3;
    
    % Append data for each trial
    if sum(qFrames) > 0
        quiescenceEvents = quiescenceEvents.append_annotation_data(expID, iBlock, qFrames, ...
                frameTimes);
    end
    if sum(iFrames) > 0
        isoMoveEvents = isoMoveEvents.append_annotation_data(expID, iBlock, iFrames, frameTimes);
    end
    if sum(gFrames) > 0
        groomEvents = groomEvents.append_annotation_data(expID, iBlock, gFrames, frameTimes);
    end
    if sum(lFrames) > 0
        locEvents = locEvents.append_annotation_data(expID, iBlock, lFrames, frameTimes);
    end
   
end%iBlock

% Reference images (using volAvgSessionData to get one set for each block)
volAvgRefImgFile = dir(fullfile(expDir, '*volAvgSessionData.mat'));
load(fullfile(expDir, volAvgRefImgFile.name), 'volAvgSessionData');
refImages = [];
for iTrial = 1:expMd.nTrials
    refImages(:, :, :, iTrial) = mean(volAvgSessionData(:, :, :, ...
            blockStartTrials(iTrial):blockEndTrials(iTrial)), 4);
end

% Also save the full experiment reference images in an array just in case I want them later
try
    fullExpRefImages = [];
    for iPlane = 1:numel(aD.refImg)
        fullExpRefImages(:, :, iPlane) = aD.refImg{iPlane};
    end
    fullExpRefImages = repmat(fullExpRefImages, 1, 1, 1, expMd.nTrials); % --> [y, x, plane, trial];
catch; end

catch ME; rethrow(ME); end

%% Generate names for ROIs in older experiments

MAX_INTENS = 700;

try
    
roiData = {};
load(fullfile(expDir, 'ROI_metadata.mat'), 'ROImetadata');
for iROI = 1:numel(ROImetadata)
    figure(iROI);clf;
    currROI = ROImetadata{iROI};
    nPlots = numel(currROI);
    subplotDims = numSubplots(nPlots);
    for iPlot = 1:nPlots
       subaxis(subplotDims(1), subplotDims(2), iPlot);
       imshow(currROI(iPlot).refImg, [0 MAX_INTENS]);
       hold on
       plot(currROI(iPlot).xi, currROI(iPlot).yi, 'linewidth', 2, 'color', 'r')
    end
end
disp(['nROIs = ', num2str(numel(ROImetadata))])

% Activate figures in reverse order so first one is on top
for iROI = numel(ROImetadata):-1:1
   figure(iROI); 
end

catch ME; rethrow(ME); end

%% Reformat ROI data to match newer experiments

roiNames = {'TypeD-R', 'ANT-R', 'Background-SMP'};

try
    
roiData = {};
load(fullfile(expDir, 'ROI_metadata.mat'), 'ROImetadata');
% allROImd = {};
% uniqueSids = unique(trialSidNums);
% for iSid = 1:numel(uniqueSids)
%     load(fullfile(expDir, ['sid_', num2str(uniqueSids(iSid)), '_ROI_metadata.mat']), 'ROImetadata');
%     allROImd{iSid} = ROImetadata;
% end

load(fullfile(expDir, 'ROI_data_avg.mat'), 'ROIDataAvg');
nROIs = size(ROIDataAvg, 3);

% Get raw average fluorescence data for each ROI and calcualte trial-based dF/F
for iTrial = 1:expMd.nTrials
    currBlockTrials = blockStartTrials(iTrial):blockEndTrials(iTrial);
    currROIData = ROIDataAvg(:, currBlockTrials, :);                 % --> [volume, shortTrial, ROI]
    currROIData = reshape(permute(currROIData, [3 2 1]), nROIs, []); % --> [ROI, volume]
    
    roiData{iTrial} = struct();
    for iROI = 1:nROIs
        currROI = ROImetadata{iROI};
        
        % Copy metadata
        roiData{iTrial}(iROI).name = roiNames{iROI};
        roiData{iTrial}(iROI).subROIs = struct();
        for iSubROI = 1:numel(currROI)
            roiData{iTrial}(iROI).subROIs(iSubROI).plane = currROI(iSubROI).plane;
            roiData{iTrial}(iROI).subROIs(iSubROI).position = [currROI(iSubROI).xi, currROI(iSubROI).yi];
        end
        
        % Copy raw avg fluorescence for current ROI
        roiData{iTrial}(iROI).rawFl = currROIData(iROI, :);
        
        % Calculate trial-based dF/F
        roiDataSorted = sort(roiData{iTrial}(iROI).rawFl);
        baselineF = median(roiDataSorted(1:round(numel(roiDataSorted) * 0.05))); % Bottom 5% as baseline
        roiData{iTrial}(iROI).dffData = (roiData{iTrial}(iROI).rawFl - baselineF) ./ baselineF;        
        
    end
end

% Get list of all unique ROI names in experiment
allROINames = {};
for iTrial = 1:numel(roiData)
    allROINames = [allROINames, {roiData{iTrial}.name}];
end
ROIList = unique(allROINames);

% Extract and concatenate all fl data throughout experiment for each ROI
rawROIData = repmat({[]}, 1, numel(ROIList)); stdDevs = [];
for iROI = 1:numel(ROIList)
    for iTrial = 1:numel(roiData)
        if sum(strcmp({roiData{iTrial}.name}, ROIList{iROI}))
            
            currROIData = roiData{iTrial}(strcmp({roiData{iTrial}.name}, ...
                    ROIList{iROI})).rawFl;
            
            FL_THRESH = 8; % Hack to prevent the inclusion of trials when the PMT shutter was closed
            stdDevs(iROI, iTrial) = std(currROIData);
            if std(currROIData) > FL_THRESH
                rawROIData{iROI} = [rawROIData{iROI}, currROIData];
            end
            
        end
    end
end

% Calculate a whole-experiment baseline F value for each ROI
for iROI = 1:numel(ROIList)
    currFlSorted = sort(rawROIData{iROI});
    currROIBaseline = median(currFlSorted(1:round(numel(currFlSorted) * 0.05))); % This is actually just the 2.5-th percentile I see?
    
    % Calculate experiment-wide dF/F for each trial
    for iTrial = 1:numel(roiData)
         if sum(strcmp({roiData{iTrial}.name}, ROIList{iROI}))
             expDffData = (roiData{iTrial}(strcmp({roiData{iTrial}.name}, ...
                    ROIList{iROI})).rawFl - currROIBaseline)./ currROIBaseline;
             roiData{iTrial}(strcmp({roiData{iTrial}.name}, ...
                    ROIList{iROI})).expDffData = expDffData;
         end
    end    
end

catch ME; rethrow(ME); end

%% ----- Save data in group analysis directory -----
try
writetable(expMd, fullfile(saveDir, [expID, '_expMetadata.csv']));
save(fullfile(saveDir, [expID, '_trialMetadata.mat']), 'trialMd');
% writetable(trialMd, fullfile(saveDir, [expID, '_trialMetadata.csv']));
writetable(daqOutputData, fullfile(saveDir, [expID, '_daqOutputData.csv']));
save(fullfile(saveDir, [expID, '_ficTracData.mat']), 'ftData'); 
save(fullfile(saveDir, [expID, '_refImages.mat']), 'refImages');
save(fullfile(saveDir, [expID, '_fullExpRefImages.mat']), 'fullExpRefImages');
% save(fullfile(saveDir, [expID, '_roiData.mat']), 'roiData');

if ~isempty(quiescenceEvents.eventData)
    quiescenceEvents.export_csv(saveDir, 'fileNamePrefix', expID);
end
if ~isempty(isoMoveEvents.eventData)
    isoMoveEvents.export_csv(saveDir, 'fileNamePrefix', expID);
end
if ~isempty(groomEvents.eventData)
    groomEvents.export_csv(saveDir, 'fileNamePrefix', expID);
end
if ~isempty(locEvents.eventData)
    locEvents.export_csv(saveDir, 'fileNamePrefix', expID);
end

if exist(fullfile(expDir, 'sidTrialCounts.csv'), 'file')
    copyfile(fullfile(expDir, 'sidTrialCounts.csv'), ...
            fullfile(saveDir, [expID, '_sidTrialCounts.csv'])); 
end
catch ME; rethrow(ME); end

%% ----- ODOR EVENT DATA -----

odorANames = {'EtOH'};
odorAConcentrations = {'neat'};
odorBNames = {''};
odorBConcentrations = {'e-1'};
flowRates = {20};
trialNums = {[]};

try 
    
if any(daqOutputData.odorA + daqOutputData.odorB)
    odorEvents = odorEvent();
    
    for iCond = 1:numel(trialNums)
        
        % Get info for current set of conditions
        currNameA = odorANames{iCond};
        currNameB = odorBNames{iCond};
        concA = odorAConcentrations{iCond};
        concB = odorBConcentrations{iCond};
        if numel(flowRates) == 1
            currFlowRate = flowRates{1};
        else
            currFlowRate = flowRates{iCond};
        end
        if isempty(trialNums{iCond})
            currTrialNums = 1:expMd.nTrials;
        else
            currTrialNums = trialNums{iCond};
        end
        mdFieldNames = {'odorName', 'concentration', 'flowRate'};
        
        for iTrial = 1:numel(currTrialNums)
            
            currTrialData = daqOutputData(daqOutputData.trialNum == currTrialNums(iTrial), :);
            daqSampleTimes = (1:size(currTrialData, 1)) / expMd.daqSampRate;
            
            % Odor A 
            if any(currTrialData.odorA)
                
                mdFieldVals = {currNameA, concA, currFlowRate}; 
                
                % Identify onset and offset samples
                onsetSamples = find(diff(currTrialData.odorA) == 1) + 1;
                if currTrialData.odorA(1) == 1
                   onsetSamples = [1; onsetSamples]; 
                end
                offsetSamples = find(diff(currTrialData.odorA) == -1) + 1;
                if currTrialData.odorA(end) == 1
                   offsetSamples = [offsetSamples; numel(currTrialData.odorA)]; 
                end
                
                % Convert to times
                onsetTimes = daqSampleTimes(onsetSamples);
                offsetTimes = daqSampleTimes(offsetSamples);
                eventTimes = {[onsetTimes', offsetTimes']};              
                
                odorEvents = odorEvents.append_data(expMd.expID{:}, currTrialNums(iTrial), ...
                        eventTimes, mdFieldNames, mdFieldVals);
            end%if
            
            % Odor B 
            if any(currTrialData.odorB)
                
                mdFieldVals = {currNameB, concB, currFlowRate}; 
                
                % Identify onset and offset samples
                onsetSamples = find(diff(currTrialData.odorB) == 1) + 1;
                if currTrialData.odorB(1) == 1
                   onsetSamples = [1, onsetSamples]; 
                end
                offsetSamples = find(diff(currTrialData.odorB) == -1) + 1;
                if currTrialData.odorB(end) == 1
                   offsetSamples = [offsetSamples, numel(currTrialData.odorB)]; 
                end
                
                % Convert to times
                onsetTimes = daqSampleTimes(onsetSamples);
                offsetTimes = daqSampleTimes(offsetSamples);
                eventTimes = {[onsetTimes', offsetTimes']};              
                
                odorEvents = odorEvents.append_data(expMd.expID{:}, currTrialNums(iTrial), ...
                        eventTimes, mdFieldNames, mdFieldVals);
            end%if

        end%iTrial
    end%iCond
    
    % Export to .csv file
    odorEvents.export_csv(saveDir, 'fileNamePrefix', expMd.expID{:});
    
end%if

catch ME; rethrow(ME); end

%% ----- OPTO STIM EVENT DATA -----

LEDpowers = {20};
dutyCycles = {100};
trialNums = {[]}; % AKA block nums

try 
    
% Opto stim
if any(daqOutputData.optoStimSmoothed)
    optoStimEvents = optoStimEvent();
    
    for iCond = 1:numel(trialNums)
        
        % Get info for current set of conditions
        currPower = LEDpowers{iCond};
        currDutyCycle = dutyCycles{iCond};
        if isempty(trialNums{iCond})
            currTrialNums = 1:expMd.nTrials;
        else
            currTrialNums = trialNums{iCond};
        end
        mdFieldNames = {'LEDpower', 'dutyCycle'};
        mdFieldValues = {currPower, currDutyCycle};
        
        for iTrial = 1:numel(currTrialNums)
            currTrialData = daqOutputData(daqOutputData.trialNum == currTrialNums(iTrial), :);
            daqSampleTimes = (1:size(currTrialData, 1)) / expMd.daqSampRate;
            optoStimSamples = currTrialData.optoStimSmoothed > 0;
            
            if any(optoStimSamples)
                                
                % Identify onset and offset samples
                onsetSamples = find(diff(optoStimSamples) == 1) + 1;
                if optoStimSamples(1) == 1
                    onsetSamples = [1; onsetSamples];
                end
                offsetSamples = find(diff(optoStimSamples) == -1) + 1;
                if optoStimSamples(end) == 1
                    offsetSamples = [offsetSamples; numel(optoStimSamples)];
                end
                
                % Convert to times
                onsetTimes = daqSampleTimes(onsetSamples);
                offsetTimes = daqSampleTimes(offsetSamples);
                eventTimes = {[onsetTimes', offsetTimes']};
                
                optoStimEvents = optoStimEvents.append_data(expMd.expID{:}, currTrialNums(iTrial), ...
                        eventTimes, mdFieldNames, mdFieldVals);
            end%if
            
        end%iTrial
        
    end%iCond
    
    % Export to .csv file
    optoStimEvents.export_csv(saveDir, 'fileNamePrefix', expMd.expID{:});
    
end

catch ME; rethrow(ME); end

%% Speaker
if any(daqOutputData.speaker)
    
end 











