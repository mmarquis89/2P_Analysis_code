%%% Closed loop analyis testing code
% 
% fictracData --> [???, odorA valve, odorB shutoff, NO valve, yaw]
% Arduino sampling rate = 4000 Hz;
%
% OdorOn voltage ranges: 
%       2018_07_07_exp_1_tid_1-2:   [1.77 3.6] V 
%       2018_07_07_exp_1_tid_3:     [1.77 5.2] V
%       2018_07_07_exp_1_tid_4:     [0 3.5] V
%       2018_07_07_exp_2:           [0 5] V

odorOnRange = [0 5];
tid = 1;

%% LOAD REGISTERED DATA FILE AND BEHAVIORAL ANNOTATION DATA


[dataFile, parentDir, ~] = uigetfile('*.mat', 'Select an imaging data session file', 'D:\Dropbox (HMS)\2P Data\Imaging Data\');


SAMP_RATE = 40;

try
    if dataFile == 0
        % Skip loading if user clicked "Cancel"
        disp('Initialization cancelled')
    else
        % Load matfile object to query session data size
        disp('Loading matfile...')
        m = matfile([parentDir, dataFile]); % Only field is 'wholeSession'
        
        % Load analysis metadata
        disp('Loading analysis metadata...')
        load(fullfile(parentDir, 'analysisMetadata.mat')) % --> 'analysisMetadata' info struct
        
        % Load reference images file
        disp('Loading reference images...')
%         [refImgFile, refImgPath, ~] = uigetfile('*.mat', 'Select a reference image file', parentDir);
%         refImgPath = parentDir; refImgFile = ['sid_', num2str(analysisMetadata.sid), '_refImages.mat'];
        refImgPath = parentDir; refImgFile = ['refImages_Reg.mat'];
        load(fullfile(refImgPath, refImgFile)) % --> 'refImages', 'channelNum'
        clear refImgPath refImgFile
        
        % Load PCA data
        disp('Loading PCA data...')
        if exist(fullfile(parentDir, ['PCA_data_', dataFile ]), 'file')
            load(fullfile(parentDir, ['PCA_data_', dataFile ])) % --> 'explained', 'pcaData' ([pc, plane], [y, x, pc, plane]
        end
        
        % Load annotation type data
        disp('Loading annotation type data...')
        load(fullfile(parentDir, 'annotationTypes.mat')) % --> 'annotationTypes', 'annotationTypeSummary'
        
        % Load FicTrac data
        disp('Loading FicTrac data...')
        ftDir = fullfile('D:\Dropbox (HMS)\2P Data\Behavior Vids\', analysisMetadata.expDate, '_Movies\FicTracData');
        if isdir(ftDir)
            ftData = load_fictrac_data(analysisMetadata, 'Sid', analysisMetadata.sid, 'ParentDir', ftDir);
        else
            ftData = load_fictrac_data(analysisMetadata, 'Sid', analysisMetadata.sid);
        end
        if ~isfield(analysisMetadata, 'ftData')
            analysisMetadata.ftData = ftData;
            save(fullfile(parentDir, 'analysisMetadata.mat'), 'analysisMetadata', '-v7.3');
        end
        
        % Load volume-averaged raw session data
        volAvgDataFile = fullfile(parentDir, ['sid_', num2str(analysisMetadata.sid), '_volAvgSessionData.mat']);
        if exist(volAvgDataFile)
            load(volAvgDataFile) % "volAvgSessionData"
        end
        
        % Reload analysis metadata in case .goodTrials has been updated
        load(fullfile(parentDir, 'analysisMetadata.mat')) % --> 'analysisMetadata' info struct
        analysisMetadata.refImg = refImages;
        disp('All data loaded')
        
        % Omit any trials in which FicTrac reset from analysis
        analysisMetadata.goodTrials(logical(ftData.resets)) = 0;
        
        % ------- Copy variables for convenience -------
        sessionSize = size(m, 'wholeSession');
        expDate = analysisMetadata.expDate;
        sid = analysisMetadata.sid;
        nPlanes = analysisMetadata.nPlanes;
        nVolumes = analysisMetadata.nVolumes;
        refImg = analysisMetadata.refImg;
        if ~isempty(analysisMetadata.nFrames)
            nFrames = analysisMetadata.nFrames;
        else
            nFrames = nVolumes;
        end
        nTrials = analysisMetadata.nTrials;
        nGoodTrials = sum(analysisMetadata.goodTrials);
        stimTypes = analysisMetadata.stimTypes;
        stimOnsetTimes = analysisMetadata.stimOnsetTimes;
        stimDurs = analysisMetadata.stimDurs;
        trialDuration = analysisMetadata.trialDuration;
        volumeRate = analysisMetadata.volumeRate;
        volFrames = analysisMetadata.volFrames;
        goodTrials = analysisMetadata.goodTrials;
        stimSepTrials = analysisMetadata.stimSepTrials;
        behaviorAnnotArr = annotationTypes{contains(annotationTypeSummary.AnnotationType, 'move')}.frameAnnotArr;
        
        % Create hardcoded parameters
        FRAME_RATE = 25; % This is the frame rate for the behavior video
        MAX_INTENSITY = analysisMetadata.MAX_INTENSITY;
        if isempty(nFrames)
            nFrames = sum(trialDuration) * FRAME_RATE;
        end
        volTimes = (1:nVolumes)' ./ volumeRate;
        frameTimes = (1:nFrames)' ./ FRAME_RATE;
        nSamps = size(analysisMetadata.daqFtData, 2);
        sampTimes = (1:nSamps)' ./ SAMP_RATE;
        volSamps = [];
        for iVol = 1:nVolumes
            [~, volSamps(iVol)] = min(abs(sampTimes - volTimes(iVol)));
        end
        % Create directory for saving analysis files if necessary
        if ~isdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'])
            mkdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'])
        end
        
        parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)];
        metaDataFileName = 'ROI_metadata.mat';
        dffDataFileName = 'ROI_Data_Avg.mat';
        
        % Load metadata
        load(fullfile(parentDir, metaDataFileName)); % --> ROImetadata(.mask, .xi, .yi, .plane, .color, .refImg)
        analysisMetadata.ROImetadata = ROImetadata;
        analysisMetadata.nROIs = numel(ROImetadata); nROIs = analysisMetadata.nROIs;
        
        % Load imaging and dF/F data
        load(fullfile(parentDir, dffDataFileName)); % --> 'ROIDataAvg', 'ROIDffAvg', 'ROIDataBaseSub' ([volume, trial, ROI])
        disp('ROI data loaded')
        
        % Plot ROIs
        plot_ROIs(ROImetadata);
        
        % Plot mean value in each ROI for each trial
        volAvgROIData = squeeze(mean(ROIDataAvg, 1)); % --> [trial, ROI]
        figure(nROIs + 1);clf;hold on
        legendStr = [];
        for iROI = 1:nROIs
            currROIData = volAvgROIData(:, iROI);
            zeroedData = currROIData - min(currROIData(:));
            legendStr{iROI} = num2str(iROI);
            plot(zeroedData);
        end
        legend(legendStr);
        
        
    end%if
    
catch foldME; rethrow(foldME); end

% Initial processing steps
daqFtData = analysisMetadata.daqFtData; % --> [var, frame, block]
daqFtData = daqFtData(:,:,tid);
stimOnData = squeeze(daqFtData(5, :, :)); % --> [frame, block]
odorOnRangeRad = odorOnRange * (1/max(stimOnData(:))) * 2*pi;

% radYaw = squeeze(daqFtData(3, :, :)) .* (1/max(stimOnData(:))) * 2*pi;
xData = squeeze(daqFtData(1, :, :)) .* (1/max(stimOnData(:))) .* 4.5;
yData = squeeze(daqFtData(2, :, :)) .* (1/max(stimOnData(:))) .* 4.5;
yawData = squeeze(daqFtData(3, :, :)) .* (1/max(stimOnData(:))) * 2*pi;


odorOn = logical(round(squeeze(daqFtData(5, :, :))));
odorOn(1) = 0;
odorOn(end) = 0;

% Calculate speed data
xSpeed = smooth([0, diff(xData)], 1) * SAMP_RATE;
ySpeed = smooth([0, diff(yData)], 1) * SAMP_RATE;
yawSpeed = rad2deg(smooth([0, diff(yawData)], 1) * SAMP_RATE);

% Remove erroneously high speeds due to unwrapping
xSpeed(abs(xSpeed) > (2 * std(xSpeed))) = nan;
ySpeed(abs(ySpeed) > (2 * std(ySpeed))) = nan;
yawSpeed(abs(yawSpeed) > (2 * std(yawSpeed))) = nan;

% Smooth speed data
smXSpeed = repeat_smooth(xSpeed, 7, 1, 2);
smYSpeed = repeat_smooth(ySpeed, 7, 1, 2);
smYawSpeed = repeat_smooth(yawSpeed, 7, 1, 2);

% Get onset/offset/event times
odorOnStr = regexprep(num2str(odorOn), ' ', '');
[onsetInds, offsetInds] = regexp(odorOnStr, '01+0');
odorOnsets = zeros(size(odorOn));
odorOffsets = zeros(size(odorOn));
odorOnsets(onsetInds + 1) = 1;
odorOffsets(offsetInds) = 1;

% Make event list of odor presentation
odorEventList = create_event_list(odorOnsets, odorOffsets);
eventDurs = odorEventList(:,2) - odorEventList(:,1);
eventDursSec = eventDurs / SAMP_RATE;
ieiList = odorEventList(2:end, 1) - odorEventList(1:end - 1, 2);
ieiList = [odorEventList(1,1); ieiList];


%% 

% Plot X (forward)
f = figure(1); clf; 
subaxis(3,1,1); hold on; 
plot(smXSpeed, '-.');
ax = gca();
ax.XTick = 1:SAMP_RATE:numel(odorOn);
ax.XTickLabel = (1:1:ceil((numel(odorOn) / SAMP_RATE)));

% Shade around odor stims
for iEvent = 1:size(odorEventList, 1)
    startSamp = odorEventList(iEvent, 1);
    endSamp = odorEventList(iEvent, 2);
    xData = [startSamp, startSamp, endSamp, endSamp];
    yData = [min(smXSpeed(:)), max(smXSpeed(:)), max(smXSpeed(:)), min(smXSpeed(:))];
    fill(xData, yData, rgb('red'), 'facealpha', 0.2, 'edgealpha', 0)
end

% Plot Y (lateral)
subaxis(3,1,2); hold on; 
plot(smYSpeed, '-.');
ax = gca();
ax.XTick = 1:SAMP_RATE:numel(odorOn);
ax.XTickLabel = (1:1:ceil((numel(odorOn) / SAMP_RATE)));

% Shade around odor stims
for iEvent = 1:size(odorEventList, 1)
    startSamp = odorEventList(iEvent, 1);
    endSamp = odorEventList(iEvent, 2);
    xData = [startSamp, startSamp, endSamp, endSamp];
    yData = [min(smYSpeed(:)), max(smYSpeed(:)), max(smYSpeed(:)), min(smYSpeed(:))];
    fill(xData, yData, rgb('red'), 'facealpha', 0.2, 'edgealpha', 0)
end

% Plot yaw data
subaxis(3,1,3); hold on; 
plot(smYawSpeed, '-.');
ax = gca();
ax.XTick = 1:SAMP_RATE:numel(odorOn);
ax.XTickLabel = (1:1:ceil((numel(odorOn) / SAMP_RATE)));

% Shade around odor stims
for iEvent = 1:size(odorEventList, 1)
    startSamp = odorEventList(iEvent, 1);
    endSamp = odorEventList(iEvent, 2);
    xData = [startSamp, startSamp, endSamp, endSamp];
    yData = [min(smYawSpeed(:)), max(smYawSpeed(:)), max(smYawSpeed(:)), min(smYawSpeed(:))];
    fill(xData, yData, rgb('red'), 'facealpha', 0.2, 'edgealpha', 0)
end


%% PLOT GCAMP RESPONSE TO EACH EVENT

% Get ROI data for each event
baselineDur = 1; respDur = 3;
currTrialROI = squeeze(ROIDataAvg(:,tid, :)); % --> [volume, ROI]
eventFlData = [];
for iEvent = 1:size(odorEventList, 1)
    if ieiList(iEvent) > baselineDur * SAMP_RATE && eventDurs(iEvent) > respDur * SAMP_RATE
        [~, startVol] = min(abs(volSamps - odorEventList(iEvent, 1)));
        [~, endVol] = min(abs(volSamps - odorEventList(iEvent, 2)));
        eventFlData{iEvent} = currTrialROI(floor(startVol-(baselineDur * volumeRate)):floor(startVol + (respDur * volumeRate)), :);
    end
end

odorEventFlArr = [];
for iEvent = 1:numel(eventFlData)
   if ~isempty(eventFlData{iEvent})
      odorEventFlArr = cat(3, odorEventFlArr, eventFlData{iEvent}); % --> [vol, ROI, event]
   end
end
odorEventFlArr = permute(odorEventFlArr, [1 3 2]); % --> [vol, event, ROI]
for iROI = 1:size(odorEventFlArr, 3)
    
    f = figure(iROI); clf;
    plot_ROI_data(gca, odorEventFlArr(:,:,iROI))
    
end

%% PLOT FICTRAC RESPONSE AROUND EACH EVENT

% Get FicTrac data for each event
baselineDur = 1; respDur = 3;
eventFtData = [];
for iEvent = 1:size(odorEventList, 1)
    if iei(iEvent) > baselineDur * SAMP_RATE && eventDurs(iEvent) > respDur * SAMP_RATE
        startSamp = odorEventList(iEvent, 1) - (baselineDur * SAMP_RATE);
        endSamp = odorEventList(iEvent, 1) + (respDur * SAMP_RATE);
        eventFtData(:, iEvent, 1) = smXSpeed(startSamp:endSamp);
        eventFtData(:, iEvent, 2) = smYawSpeed(startSamp:endSamp); % --> [samp, event, var]
    end
end

figure(1);clf

