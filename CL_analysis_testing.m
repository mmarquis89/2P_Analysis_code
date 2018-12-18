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
        
        % Create directory for saving analysis files if necessary
        if ~isdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'])
            mkdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'])
        end
        
    end%if
catch foldME; rethrow(foldME); end

FRAME_RATE = 40;

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
xSpeed = smooth([0, diff(xData)], 1) * FRAME_RATE;
ySpeed = smooth([0, diff(yData)], 1) * FRAME_RATE;
yawSpeed = rad2deg(smooth([0, diff(yawData)], 1) * FRAME_RATE);

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
eventDursSec = eventDurs / FRAME_RATE;



%%

% Plot X
f = figure(1); clf; 
subaxis(3,1,1); hold on; 
plot(smXSpeed, '-.');
ax = gca();
ax.XTick = 1:FRAME_RATE:numel(odorOn);
ax.XTickLabel = (1:1:ceil((numel(odorOn) / FRAME_RATE)));

% Shade around odor stims
for iEvent = 1:size(odorEventList, 1)
    startSamp = odorEventList(iEvent, 1);
    endSamp = odorEventList(iEvent, 2);
    xData = [startSamp, startSamp, endSamp, endSamp];
    yData = [min(smXSpeed(:)), max(smXSpeed(:)), max(smXSpeed(:)), min(smXSpeed(:))];
    fill(xData, yData, rgb('red'), 'facealpha', 0.2, 'edgealpha', 0)
end

% Plot X
subaxis(3,1,2); hold on; 
plot(smYSpeed, '-.');
ax = gca();
ax.XTick = 1:FRAME_RATE:numel(odorOn);
ax.XTickLabel = (1:1:ceil((numel(odorOn) / FRAME_RATE)));

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
ax.XTick = 1:FRAME_RATE:numel(odorOn);
ax.XTickLabel = (1:1:ceil((numel(odorOn) / FRAME_RATE)));

% Shade around odor stims
for iEvent = 1:size(odorEventList, 1)
    startSamp = odorEventList(iEvent, 1);
    endSamp = odorEventList(iEvent, 2);
    xData = [startSamp, startSamp, endSamp, endSamp];
    yData = [min(smYawSpeed(:)), max(smYawSpeed(:)), max(smYawSpeed(:)), min(smYawSpeed(:))];
    fill(xData, yData, rgb('red'), 'facealpha', 0.2, 'edgealpha', 0)
end












fwSmooth = smooth(ftVidData.fwSpeed(:,tid), 3);
fwNorm = fwSmooth / max(fwSmooth);
plot(fwNorm)

ax = gca();
ax.XTick = 1:FRAME_RATE:numel(odorOn);
ax.XTickLabel = (1:1:ceil((numel(odorOn) / FRAME_RATE)));

% Shade around odor stims
uwYaw = unwrap(radYaw);
for iEvent = 1:size(odorEventList, 1)
    startSamp = odorEventList(iEvent, 1);
    endSamp = odorEventList(iEvent, 2);
    xData = [startSamp, startSamp, endSamp, endSamp];
    yData = [0 1 1 0];
    fill(xData, yData, rgb('red'), 'facealpha', 0.5, 'edgealpha', 0)
end

% % Plot yaw around every event
% for iEvent = 1:size(odorEventList, 1)
% %    figure(iEvent); clf; hold on;
%    
%    
%     
% end


