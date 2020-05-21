saveDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
expList = load_expList('groupname', 'singleTrialAcq_preFicTrac');
expList = expList(3:end, :);

expNum =  28;

% ----- AUTOMATIC PROCESSING AND CONVERSION -----
try
    
% Convert expID to old-format expDir name
currExpID = expList.expID{expNum};
disp(currExpID);
expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)];
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

FRAME_RATE = 25; % Behavior video frame rate was a constant;

% ----- CREATE EXP METADATA TABLE -----
expID = currExpID;
if ~strcmp(expID(end - 1), '-')
   expID = [expID, '-1']; % Add exp num of 1 if there wasn't one in the expDate
end
expMd = table({expID}, 'VariableNames', {'expID'});
expMd.expName = {expName};
expMd.daqSampRate = nan;
expMd.panelsDisplayRate = 50;
expMd.volumeRate = 6.44;    % This is true for all these early experiments
expMd.nPlanes = 12;         % Same deal for this
annotFile = dir(fullfile(expDir, 'sid*Annotations.mat'));
if numel(annotFile) > 1
    error('Multiple annotation files found');
end
load(fullfile(expDir, annotFile(1).name), 'trialAnnotations'); 
expMd.nTrials = numel(trialAnnotations);

% ----- CREATE TRIAL METADATA TABLE -----
imgDataFiles = dir(fullfile(expParentDir, expDirName, 'cdata*.mat'));
m = matfile(fullfile(expParentDir, expDirName, imgDataFiles(1).name));
sz = size(m, 'imgData');
trialMd = [];
for iTrial = 1:expMd.nTrials
    newRow = table({expID}, 'VariableNames', {'expID'});
    newRow.trialNum = iTrial;
    newRow.trialDuration = size(trialAnnotations{1}, 1) / FRAME_RATE;
    newRow.nVolumes = sz(4);
    newRow.nDaqSamples = nan;
    newRow.nPanelsFrames = expMd.panelsDisplayRate * newRow.trialDuration;
    newRow.usingOptoStim = 0;
    newRow.usingPanels = 0;
    newRow.using2P = 1;
    trialMd = [trialMd; newRow];
end

% ----- BEHAVIOR ANNOTATION DATA -----

% Load annotation data file 
load(fullfile(expDir, annotFile(1).name), 'trialAnnotations'); 

quiescenceEvents = behaviorEvent('Quiescence');
isoMoveEvents = behaviorEvent('IsolatedMovement');
groomEvents = behaviorEvent('Grooming');
locEvents = behaviorEvent('Locomotion');
ballStopEvents = behaviorEvent('BallStop');
for iTrial = 1:expMd.nTrials

        currAnnotData = trialAnnotations{iTrial};
        frameTimes =  currAnnotData.frameTime;
        
        % Separate according to behavior type
        qFrames = currAnnotData.actionNums == 0;
        iFrames = currAnnotData.actionNums == 4;
        gFrames = currAnnotData.actionNums == 3;
        lFrames = currAnnotData.actionNums == 2;
        if any(strcmp(fieldnames(currAnnotData), 'ballStopNums'))
            bFrames = logical(currAnnotData.ballStopNums);
            bFrames(1) = 0;
        else
           bFrames = []; 
        end
        
        % Append data for each trial
        if sum(qFrames) > 0
            quiescenceEvents = quiescenceEvents.append_annotation_data(expID, iTrial, qFrames, ...
                    frameTimes);
        end
        if sum(iFrames) > 0
            isoMoveEvents = isoMoveEvents.append_annotation_data(expID, iTrial, iFrames, ...
                    frameTimes);
        end
        if sum(gFrames) > 0
            groomEvents = groomEvents.append_annotation_data(expID, iTrial, gFrames, frameTimes);
        end
        if sum(lFrames) > 0
            locEvents = locEvents.append_annotation_data(expID, iTrial, lFrames, frameTimes);
        end
        if sum(bFrames) > 0
            ballStopEvents = ballStopEvents.append_annotation_data(expID, iTrial, bFrames, ...
                    frameTimes);
        end
end%iTrial


% ----- REFERENCE IMAGES -----

% Just using the reference images from the first block as a representative sample
refImageFiles = dir(fullfile(expParentDir, expDirName, 'refImages_reg*.mat'));
load(fullfile(expParentDir, expDirName, refImageFiles(1).name), 'refImages');
fullExpRefImages = refImages;

catch ME; rethrow(ME); end


%% ----- Save data in group analysis directory -----

try
writetable(expMd, fullfile(saveDir, [expID, '_expMetadata.csv']));
writetable(trialMd, fullfile(saveDir, [expID, '_trialMetadata.csv']));
% save(fullfile(saveDir, [expID, '_refImages.mat']), 'refImages');
save(fullfile(saveDir, [expID, '_fullExpRefImages.mat']), 'fullExpRefImages');

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
if ~isempty(ballStopEvents.eventData)
    ballStopEvents.export_csv(saveDir, 'fileNamePrefix', expID);
end

if exist(fullfile(expDir, 'sidTrialCounts.csv'), 'file')
    copyfile(fullfile(expDir, 'sidTrialCounts.csv'), ...
            fullfile(saveDir, [expID, '_sidTrialCounts.csv'])); 
end
catch ME; rethrow(ME); end


%% ----- ODOR/LASER EVENT DATA -----

odorANames = {'ACV'};
odorAConcentrations = {'e-2'};
odorBNames = {'NULL'};
odorBConcentrations = {'neat'};
flowRates = {30};

try 
odorEvents = odorEvent();
laserEvents = stimEvent('IR-Laser');
mdFiles = dir(fullfile(expDir, 'metadata*.mat'));
if numel(mdFiles) ~= expMd.nTrials
   error('metadata nTrials mismatch'); 
end
for iFile = 1:numel(mdFiles)
    currFile = fullfile(expDir, mdFiles(iFile).name);
    disp(mdFiles(iFile).name);
    load(currFile, 'metaData');
    if any(strcmp(fieldnames(metaData), 'odorStimTiming'))
        odorStimTiming = metaData.odorStimTiming;
        trialDuration = metaData.trialDuration;
        if strcmp(metaData.stimType, 'OdorA') 
            preStim = trialDuration(1);
            stimDur = odorStimTiming(1);
            ISI = odorStimTiming(3);
            odorEvents = odorEvents.append_shorthand(currExpID, iFile, [preStim, stimDur, ISI], ...
                    sum(trialDuration(1:2)), odorANames{:}, odorAConcentrations{:}, flowRates{:});
                
        elseif strcmp(metaData.stimType, 'OdorB') || strcmp(metaData.stimType, 'OdorBPulse')
            preStim = trialDuration(1);
            stimDur = odorStimTiming(1);
            ISI = odorStimTiming(3);
            odorEvents = odorEvents.append_shorthand(currExpID, iFile, [preStim, stimDur, ISI], ...
                sum(trialDuration(1:2)), odorBNames{:}, odorBConcentrations{:}, flowRates{:});           
        end
    elseif strcmp(metaData.stimType, 'CenterWind')
        preStim = trialDuration(1);
        stimDur = trialDuration(2);
        ISI = trialDuration(3);
        odorEvents = odorEvents.append_shorthand(currExpID, iFile, [preStim, stimDur, ISI], ...
                sum(trialDuration), 'CleanAir', 'neat', flowRates{:});           
            
    elseif numel(metaData.trialDuration) == 1 
        preStim = str2double(regexp(metaData.stimType, '(?<=Onset-).*(?=-Dur)', 'match', 'once'));
        stimDur = str2double(regexp(metaData.stimType, '(?<=Dur-).*', 'match', 'once'));
        ISI = metaData.trialDuration;
        
        % Process any odor/carrier stream stop events
        if contains(metaData.stimType, 'OdorA')
            stimName = odorANames{:};
            stimConc = odorAConcentrations{:};
        elseif contains(metaData.stimType, 'OdorB')
            stimName = odorBNames{:};
            stimConc = odorBConcentrations{:};
        elseif contains(metaData.stimType, 'CarrierStreamStop')
            stimName = 'CarrierStreamStop';
            stimConc = 'neat';
        end
        
        % Process any IR laser events
        if contains(metaData.stimType, 'Laser')
            laserEvents = laserEvents.append_shorthand(currExpID, iFile, [preStim, stimDur, ISI], ...
                    metaData.trialDuration);
        else
            odorEvents = odorEvents.append_shorthand(currExpID, iFile, [preStim, stimDur, ISI], ...
                    metaData.trialDuration, stimName, stimConc, flowRates{:});
        end
        
    end
    
end 

% Save event data
if ~isempty(odorEvents.eventData)
    odorEvents.export_csv(saveDir, 'fileNamePrefix', expMd.expID{:});
end
if ~isempty(laserEvents.eventData)
    laserEvents.export_csv(saveDir, 'fileNamePrefix', expMd.expID{:});
end


catch ME; rethrow(ME); end

disp('Odor event data saved');

%% ----- SPEAKER EVENT DATA -----

toneFreqs = {200};
ampSettings = {-30};
trialNums = {[]};

try 
    
if any(daqOutputData.speaker)
    soundEvents = soundStimEvent();
    
    for iCond = 1:numel(trialNums)
        
        % Get info for current set of conditions
        if numel(toneFreqs) == 1
            currFreq = toneFreqs{1};
        else
            currFreq = toneFreqs{iCond};
        end
        currAmpSetting = ampSettings{iCond};
        if isempty(trialNums{iCond})
            currTrialNums = 1:expMd.nTrials;
        else
            currTrialNums = trialNums{iCond};
        end
        mdFieldNames = {'toneFreq', 'ampSetting'};
        
        for iTrial = 1:numel(currTrialNums)
            
            currTrialData = daqOutputData(daqOutputData.trialNum == currTrialNums(iTrial), :);
            daqSampleTimes = (1:size(currTrialData, 1)) / expMd.daqSampRate;
            
            if any(currTrialData.speaker)
                
                mdFieldVals = {currFreq, currAmpSetting}; 
                
                % Identify onset and offset samples
                onsetSamples = find(diff(logical(currTrialData.speaker)) == 1) + 1;
                if logical(currTrialData.speaker(1)) == 1
                   onsetSamples = [1; onsetSamples]; 
                end
                offsetSamples = find(diff(logical(currTrialData.speaker)) == -1) + 1;
                if logical(currTrialData.speaker(end)) == 1
                   offsetSamples = [offsetSamples; numel(currTrialData.speaker)]; 
                end
                
                % Convert to times
                onsetTimes = daqSampleTimes(onsetSamples);
                offsetTimes = daqSampleTimes(offsetSamples);
                eventTimes = {[onsetTimes', offsetTimes']};              
                
                soundEvents = soundEvents.append_data(expMd.expID{:}, currTrialNums(iTrial), ...
                        eventTimes, mdFieldNames, mdFieldVals);
            end%if
        end%iTrial
    end%iCond
    
    % Export to .csv file
    soundEvents.export_csv(saveDir, 'fileNamePrefix', expMd.expID{:});
    
end%if

catch ME; rethrow(ME); end










