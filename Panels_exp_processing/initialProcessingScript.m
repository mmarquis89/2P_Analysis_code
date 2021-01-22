startDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';

expDir = uigetdir(startDir, 'Select an experiment directory');
outputDir = fullfile(expDir, 'ProcessedData');

if ~isfolder(outputDir)
    mkdir(outputDir)
end

%% Make anatomy stack
try
create_anatomy_stack(expDir, 'FileString', 'Stack_*.tif', 'OutputFilePrefix', 'AnatomyStack', ...
        'OutputDir', outputDir);
catch
   disp('Error: anatomy stack creation failed'); 
end

%% Process raw imaging data
try
    
% Identify imaging data files
imgDataFiles = dir(fullfile(expDir, '*trial*.tif'));

medVals = {};
minVals = {};
minFrameVals = {};
volCounts = [];
for iFile = 1:numel(imgDataFiles)
    
    % Load the current file
    currFile = fullfile(expDir, imgDataFiles(iFile).name);
    disp(currFile);
    [imgData, siMetadata] = read_tif(currFile, 0); % --> [y, x, planes, volumes]
    sz = size(imgData);
    volCounts(iFile) = size(imgData, 4);
    
    % Extract the trial number from the fileName
    trialNum = get_trialNum(currFile);
    
    % Detect the number of channels, if there are >1 process only the first channel
    if numel(sz) == 5
        imgData = imgData(:, :, :, :, 1);
    end    
    
    % Reformat ScanImage metadata
    siMetadata = parse_scanimage_metadata(siMetadata);
    
    % Discard flyback frames
    nFlybackFrames = siMetadata.SI.hFastZ.numDiscardFlybackFrames;
    imgData(:, :, (end - (nFlybackFrames - 1)):end, :) = []; % --> [y, x, plane, volume]
    
    % Calculate some volume-by-volume summary statistics for comparison across trials
    minFrameVals{iFile} = sort(as_vector(min(min(imgData))));
    rsData = reshape(permute(imgData, [4 1 2 3]), size(imgData, 4), []); % --> [volume, px]
    medVals{iFile} = median(rsData, 2);
    minVals{iFile} = min(rsData, [], 2);
    
    % Save raw imaging data as a .mat file
    save(fullfile(outputDir, ['imagingData_raw_trial_', pad(num2str(trialNum), 3, 'left', '0'), ...
            '.mat']), 'imgData');
        
    % Also save reference images
    refImages = mean(imgData, 4); % --> [y, x, plane]
    save(fullfile(outputDir, ['refImages_raw_trial_', pad(num2str(trialNum), 3, 'left', '0'), ...
            '.mat']), 'refImages');
end

% Save summary statistics and list of volume counts for all trials
save(fullfile(outputDir, 'summaryStats.mat'), 'minFrameVals', 'medVals', 'minVals');
save(fullfile(outputDir, 'volCounts.mat'), 'volCounts');

catch ME; rethrow(ME); end

%% Plane-wise NoRMCorre motion correction
try
    
% Identify raw imaging data files
imgDataFiles = dir(fullfile(outputDir, 'imagingData_raw_trial*.mat'));

for iFile = 1:numel(imgDataFiles)
    
    % Load the current file
    currFile = fullfile(outputDir, imgDataFiles(iFile).name);
    disp(currFile);
    load(currFile, 'imgData');
    
    % Motion correct each plane one at a time
    imgDataReg = zeros(size(imgData));
    regTemplates = [];
    for iPlane = 1:size(imgDataReg, 3)
        
        % Isolate data from the current plane
        currData = squeeze(imgData(:,:, iPlane, :)); % --> [y, x, volume]
        
        % Clip highest 0.1% of values and smooth with 2D gaussian filter
        srt = sort(imgData(:));
        capVal = srt(numel(srt) - round(numel(srt)/1000));
        currData(currData > capVal) = capVal;
        currDataSm = zeros(size(currData));
        for iVol = 1:size(currData, 3)
            currDataSm(:, :, iVol) = imgaussfilt(currData(:, :, iVol), 0.5);
        end
        currData = currDataSm;
        
        % Set NoRMCorre options
        options_rigid = NoRMCorreSetParms('d1', size(currData, 1), 'd2', ...
                size(currData, 2), ...
                'max_shift', [25, 25], ...
                'init_batch', 100, ...
                'us_fac', 50 ...
                );
        
        % Run registration
        [planeData, ~, regTemplate, ~] = normcorre_batch(currData, options_rigid);
        
        % Add data to full volume arrays
        imgDataReg(:, :, iPlane, :) = planeData;    % --> [y, x, plane, volume]
        regTemplates(:, :, iPlane) = regTemplate;   % --> [y, x, plane]
        
    end%iPlane
    
     imgData = imgDataReg;
     
     % Save reference images
     refImages = mean(imgData, 4); % --> [y, x, plane]
     save(regexprep(currFile, 'imagingData_raw', 'refImages_reg'), 'refImages', 'regTemplates');
     
     % Convert back to int16 data type to save space
     imgData = int16(imgData);
     
     % Save registered imaging data
     disp('Saving registered data...');
     save(regexprep(currFile, 'raw', 'reg'), 'imgData');    
     disp('Saving complete');
     
end%iFile

% Scan file names in experiment directory to get the expID
imgDataFiles = dir(fullfile(expDir, ['*trial*.tif']));
expID = imgDataFiles(1).name(1:10);

% Also create and save a consolidated refImages file
refImgFiles = dir(fullfile(outputDir, 'refImages_reg_trial_*.mat'));
allRefImages = [];
for iFile = 1:numel(refImgFiles)
    load(fullfile(outputDir, refImgFiles(iFile).name), 'refImages');
    allRefImages = cat(4, allRefImages, refImages);
end
refImages = allRefImages;
save(fullfile(outputDir, [expID, '_refImages.mat']), 'refImages');

catch ME; rethrow(ME); end

%% Create metadata tables
try

% Identify imaging data files
imgDataFiles = dir(fullfile(expDir, '*trial*.tif'));

expMd = [];
trialMd = [];
panelsMd = [];
for iFile = 1:numel(imgDataFiles)
    
    % Load ScanImage metadata for the current trial
    currFile = fullfile(expDir, imgDataFiles(iFile).name);
    disp(currFile);
    [~, siMetadata] = read_tif(currFile, 0);
    
    % Extract the trial number and expID from the fileName
    trialNum = get_trialNum(imgDataFiles(iFile).name);
    
    % Reformat ScanImage metadata
    siMetadata = parse_scanimage_metadata(siMetadata);
    
    % Load Wombat metadata for the current trial
    mdFile = dir(fullfile(expDir, ['*metadata*trial_', pad(num2str(trialNum), 3, 'left', '0'), ...
            '.mat']));
    load(fullfile(expDir, mdFile.name), 'mD'); % [trialSettings][expID][expDirName][expDir][trialNum][SAMPLING_RATE]
    
    % Load DAQ output data for the current trial
    daqFile = dir(fullfile(expDir, ['*daqData*trial_', pad(num2str(trialNum), 3, 'left', '0'), ...
            '.mat']));
    load(fullfile(expDir, daqFile.name), 'outputData', 'trialData', 'columnLabels');
    
    % Flatten structure so the contents of trialSettings are at the root level
    mD = rmfield(setstructfields(mD, mD.trialSettings), 'trialSettings');
    
    % Calculate panels information if available
    if mD.usingPanels
        xPosFunc = mD.xDimPosFun.func;
        yPosFunc = mD.yDimPosFun.func;
        panelsCycleFrames = numel(mD.xDimPosFun.func);
        panelsCycleTime = panelsCycleFrames / double(mD.displayRate);
        xPosRep = repmat(xPosFunc, 1, ceil(mD.trialDuration / panelsCycleTime));
        yPosRep = repmat(yPosFunc, 1, ceil(mD.trialDuration / panelsCycleTime));
        nPanelsFrames = mD.displayRate * mD.trialDuration;
        panelsPosX = xPosRep(1:nPanelsFrames);
        panelsPosY = yPosRep(1:nPanelsFrames);
    else
        xPosFunc = [];
        yPosFunc = [];
        panelsPosX = [];
        panelsPosY = [];
        nPanelsFrames = 0;
    end
    
    % Separate panels metadata
    pMdRow = table({mD.expID}, 'VariableNames', {'expID'});
    pMdRow.trialNum = trialNum;
    pMdRow.panelsMode = {mD.panelsMode};
    pMdRow.pattern = {mD.pattern};
    pMdRow.xDimPosFunc = {xPosFunc};
    pMdRow.yDimPosFunc = {yPosFunc};
    pMdRow.panelsPosX = {panelsPosX};
    pMdRow.panelsPosY = {panelsPosY};
    
    % Separate trial metadata
    tMd = [];
    tMd.expID = mD.expID;
    tMd.trialNum = trialNum;
    tMd.trialDuration = mD.trialDuration;
    tMd.nVolumes = siMetadata.SI.hFastZ.numVolumes;
    tMd.nDaqSamples = size(trialData, 1);
    tMd.nPanelsFrames = nPanelsFrames;
    tMd.usingOptoStim = mD.usingOptoStim;
    tMd.optoStimTiming = {mD.optoStimTiming};
    tMd.usingPanels = double(mD.usingPanels);
    tMd.using2P = double(mD.using2P);
    tMd.originalTrialCount = 1;
    tMd.pmtShutoffVols = [];
    
    % Create expMetadata struct
    if iFile == 1
        expMd.expID = mD.expID;
        expMd.expName = mD.expName;
        expMd.daqSampRate = mD.SAMPLING_RATE;
        expMd.panelsDisplayRate = double(mD.displayRate);
        expMd.volumeRate = siMetadata.SI.hRoiManager.scanVolumeRate;
        expMd.nPlanes = siMetadata.SI.hStackManager.numSlices;
    elseif isempty(expMd.panelsDisplayRate) && ~isempty(mD.displayRate)
        % If panels were used in at least one trial but not the first, update expMd data
        expMd.panelsDisplayRate = mD.displayRate;
    end
    
    % Append to metadata structs
    trialMd = [trialMd, tMd];
    panelsMd = [panelsMd; pMdRow];
    
end

expMd.nTrials = numel(trialMd);

% Convert from structs to tables
expMd = struct2table(expMd, 'AsArray', 1);
trialMetadata = struct2table(trialMd);
panelsMetadata = panelsMd;

% % Save metadata files
writetable(expMd, fullfile(outputDir, [expMd.expID{:}, '_expMetadata.csv']));
save(fullfile(outputDir, [expMd.expID{:}, '_trialMetadata.mat']), 'trialMetadata');
save(fullfile(outputDir, [expMd.expID{:}, '_panelsMetadata.mat']), 'panelsMetadata');

catch ME; rethrow(ME); end

%% Process FicTrac data 
try

ftDir = fullfile(expDir, 'FicTracData');

% Identify all FicTrac output files
ftVidFiles = dir(fullfile(ftDir, '*trial*.mp4'));
ftDataFiles = dir(fullfile(ftDir, '*trial*.dat'));
ftFrameLogFiles = dir(fullfile(ftDir, '*vidLogFrames*trial*.txt'));

% Display warning if the number of files is non consistent across types
if numel(unique([numel(ftVidFiles), numel(ftDataFiles), numel(ftFrameLogFiles)])) ~= 1
    errordlg('Warning: inconsistent file counts');
end

% Process each trial
currFtData = []; ftData = [];
for iFile = 1:numel(ftVidFiles)
    
    % Get current trial number
    trialNumStr = regexp(ftVidFiles(iFile).name, '(?<=trial_)...', 'match', 'once');
    trialNum = str2double(trialNumStr);
    disp(['Processing FicTrac data for trial #', num2str(trialNum)]);
    
    % Load output video and extract the median luminance from each frame
    currFile = fullfile(ftDir, ftVidFiles(iFile).name);
    try
        currVid = VideoReader(currFile);
    catch
       disp(['Video read error! Skipping trial #', num2str(trialNum)])
       continue 
    end
    vidData = []; medLum = [];
    frameCount = 0;
    while hasFrame(currVid)
        frameCount = frameCount + 1;
        if ~mod(frameCount, 1000)
            disp(['Reading frame ', num2str(frameCount), '...'])
        end
        currFrame = uint16(readFrame(currVid));
        medLum(frameCount) = median(as_vector(currFrame(:,:,1)));
        vidData(:,:, frameCount) = currFrame(:,:,1);
    end

    % Determine which video frames mark the beginning and end of the trial
    baseSubLum = medLum - median(medLum);
    lumThresh = 5 * std(baseSubLum) + mean(baseSubLum); % Set threshold at 4 SDs above average
    startVidFrame = find(baseSubLum > lumThresh, 1); % Index of first frame with >90% of max luminance
    endVidFrame = find(baseSubLum > lumThresh, 1, 'last');
    
    % Load vid frame log and identify first and last FicTrac data frames from the trial period
    frameLog = csvread(fullfile(ftDir, ftFrameLogFiles(iFile).name));
    startFtFrame = frameLog(startVidFrame) + 1; % Adding one because original is zero-indexed
    endFtFrame = frameLog(endVidFrame) + 1;
    
    % Load main FicTrac data file
    rawFtData = csvread(fullfile(ftDir, ftDataFiles(iFile).name));
    
    % Deal with dropped frames by repeating last recorded value
    for iFrame = 2:(size(rawFtData, 1) - 1)
        if rawFtData(iFrame, 1) ~= (rawFtData(iFrame - 1, 1) + 1)
            
            % Integrated XY position
            rawFtData(iFrame, 15:16) = rawFtData(iFrame + 1, 15:16); 
            rawFtData(iFrame:end, 15:16) = rawFtData(iFrame:end, 15:16) + ...
                    rawFtData(iFrame - 1, 15:16);
                
            % Heading direction
            rawFtData(iFrame, 17) = rawFtData(iFrame + 1, 17);
            rawFtData(iFrame:end, 17) = mod(rawFtData(iFrame:end, 17) + ...
                    rawFtData(iFrame - 1, 17), 2 * pi);
            
            % Movement speed
            rawFtData(iFrame, 19) = rawFtData(iFrame + 1, 19);
            
        end
    end
    
    % Pull out data from within the trial period
    ftTrialFrames = rawFtData(:, 1) >= startFtFrame & rawFtData(:, 1) <= endFtFrame;
    currFtData = rawFtData(ftTrialFrames, :);%rawFtData(startFtFrame:endFtFrame, :);
    currFtData(:, 1) = currFtData(:,1) - startFtFrame; % Align FrameCount to trial start
    
    % Calculate frame times for the trial video
    trialVidFrames = false(size(frameLog));
    trialVidFrames(startVidFrame:endVidFrame) = 1;
    rawFtFrameTimes = cumsum(rawFtData(:, 24)) ./ 1e9;
    rawVidFrameTimes = rawFtFrameTimes(frameLog + 1);
    
    % Save processed data files along with luminance values and frame log
    ftData(iFile).trialNum = trialNum;
    ftData(iFile).rawData = rawFtData;
    ftData(iFile).trialData = currFtData;
    ftData(iFile).frameLog = frameLog;
    ftData(iFile).medLum = medLum;
    ftData(iFile).startVidFrame = startVidFrame;
    ftData(iFile).endVidFrame = endVidFrame;
    ftData(iFile).rawVidFrameTimes = rawVidFrameTimes;
    ftData(iFile).lumThresh = lumThresh;
    
    % Write video data from within the trial period to a new file    
    vidData = vidData(:, :, trialVidFrames);
    trialVid = VideoWriter(fullfile(outputDir, ['FicTrac_video_trial_', trialNumStr]), 'MPEG-4');
    trialVid.FrameRate = round(median(1 ./ diff(rawVidFrameTimes)));
    open(trialVid)
    for iFrame = 1:size(vidData, 3)
        if ~mod(iFrame, 1000)
            disp(['Writing frame ', num2str(iFrame), '...'])
        end
        writeVideo(trialVid, uint8(vidData(:,:, iFrame)));
    end
    close(trialVid)

end%iFile
disp('Video processing complete');

% Save initial processing in raw data folder
save(fullfile(expDir, 'rawFicTracData.mat'), 'ftData');

catch ME; rethrow(ME); end

%% Define optic flow ROI around the fly and/or ball
select_video_ROIs(outputDir);

%% Calculate optic flow in FLY ROI in FicTrac videos
try 
    
% Define an ROI around the fly
roiDataFile = fullfile(outputDir, 'Behavior_Vid_ROI_Data.mat');
if exist(roiDataFile, 'file')
    load(roiDataFile, 'roiData');
    vidRoi = roiData;
else
    vidRoi = select_video_ROIs(outputDir);
end

% Extract mean flow within ROI for each FicTrac vid
ftVids = dir(fullfile(outputDir, 'FicTrac*.mp4'));
meanFlowMags = {};
validTrialNums = [];
parfor iTrial = 1:numel(ftVids)
    disp(ftVids(iTrial).name);
    tic
    trialNum = get_trialNum(ftVids(iTrial).name);
    meanFlowMags{iTrial} = optic_flow_calc(fullfile(outputDir, ftVids(iTrial).name), 'roiMask', ...
            vidRoi);
    validTrialNums(iTrial) = trialNum;
    disp(['Flow calculation completed in ', num2str(toc, '%.1f'), ' sec']);
end

% Load trialMetadata file to get a full trial list 
trialMdFile = dir(fullfile(outputDir, '*trialMetadata.mat'));
load(fullfile(outputDir, trialMdFile.name), 'trialMetadata');
trialList = [trialMetadata.trialNum];

% Adjust indexing if data is missing for any trials
newFlowMags = {};
if numel(validTrialNums) ~= numel(trialList)
    for iTrial = 1:numel(trialList)
        if ismember(iTrial, validTrialNums)
            newFlowMags{iTrial} = meanFlowMags{validTrialNums == iTrial};
        else
            newFlowMags{iTrial} = [];
        end
    end
    meanFlowMags = newFlowMags;
end

% Save flow data
save(fullfile(outputDir, 'flowMags.mat'), 'meanFlowMags')

catch ME; rethrow(ME); end

%% Calculate optic flow in BALL ROI in FicTrac videos
try 
    
% Define an ROI around the ball
roiDataFile = fullfile(outputDir, 'Behavior_Vid_ROI_Data_Ball.mat');
if exist(roiDataFile, 'file')
    load(roiDataFile, 'roiData');
    vidRoi = roiData;
else
    vidRoi = select_video_ROIs(outputDir);
end

% Extract mean flow within ROI for each FicTrac vid
ftVids = dir(fullfile(outputDir, 'FicTrac*.mp4'));
meanFlowMags = {};
validTrialNums = [];
parfor iTrial = 1:numel(ftVids)
    disp(ftVids(iTrial).name);
    tic
    trialNum = get_trialNum(ftVids(iTrial).name);
    meanFlowMags{iTrial} = optic_flow_calc(fullfile(outputDir, ftVids(iTrial).name), 'roiMask', ...
            vidRoi);
    validTrialNums(iTrial) = trialNum;
    disp(['Flow calculation completed in ', num2str(toc, '%.1f'), ' sec']);
end

% Load trialMetadata file to get a full trial list 
trialMdFile = dir(fullfile(outputDir, '*trialMetadata.mat'));
load(fullfile(outputDir, trialMdFile.name), 'trialMetadata');
trialList = [trialMetadata.trialNum];

% Adjust indexing if data is missing for any trials
newFlowMags = {};
if numel(validTrialNums) ~= numel(trialList)
    for iTrial = 1:numel(trialList)
        if ismember(iTrial, validTrialNums)
            newFlowMags{iTrial} = meanFlowMags{validTrialNums == iTrial};
        else
            newFlowMags{iTrial} = [];
        end
    end
    meanFlowMags = newFlowMags;
end

% Save flow data
save(fullfile(outputDir, 'flowMags_ballROI.mat'), 'meanFlowMags')

catch ME; rethrow(ME); end

%% Additional FicTrac processing
try

% Scan file names in experiment directory to get the expID
imgDataFiles = dir(fullfile(expDir, ['*trial*.tif']));
expID = imgDataFiles(1).name(1:10);

% Load FT vid optic flow data
clear meanFlowMagsBallROI;
ballRoiFlowFile = fullfile(outputDir, 'flowMags_ballROI.mat');
if exist(ballRoiFlowFile, 'file')
    load(ballRoiFlowFile, 'meanFlowMags')
    meanFlowMagsBallROI = meanFlowMags;
end
load(fullfile(outputDir, 'flowMags.mat'), 'meanFlowMags')

% Load raw FicTrac data
load(fullfile(expDir, 'rawFicTracData.mat'), 'ftData');

% Create table of final analysis data
rawFt = ftData;
ftData = [];
for iTrial = 1:numel(rawFt)
    currFt = rawFt(iTrial).trialData;
    
    if ~isempty(currFt)
        newRow = table({expID}, rawFt(iTrial).trialNum, 'VariableNames', {'expID', 'trialNum'});
        
        % Get frame times relative to start of trial (column 24 is nanosec since previous frame)
        ftFrameTimes = cumsum(currFt(:, 24)) ./ 1e9;
        
        % Copy important variables, converting units as needed
        IFI = [ftFrameTimes(1); diff(ftFrameTimes)];
        newRow.intX = {currFt(:, 15) * 4.5};               % mm
        newRow.intY = {currFt(:, 16) * 4.5};               % mm
        newRow.intHD = {currFt(:, 17)};                    % radians
        newRow.moveSpeed = {(currFt(:, 19) * 4.5) ./ IFI}; % mm/sec
        newRow.intFwMove = {currFt(:, 20) * 4.5};          % mm
        newRow.intSideMove = {currFt(:, 21) * 4.5};        % mm
        
        % Calculate derived FicTrac variables
        newRow.yawSpeed = {[0; diff(smoothdata(unwrap(newRow.intHD{:}, [], 1), 1, 'gaussian', 7), 1)] ...
                ./ IFI}; % radians/sec
        newRow.fwSpeed = {[0; diff(smoothdata(newRow.intFwMove{:}, 1, 'gaussian', 7), 1)] ./ IFI};     % mm/sec
        newRow.sideSpeed = {[0; diff(smoothdata(newRow.intSideMove{:}, 1, 'gaussian', 7), 1)] ...
                ./ IFI}; % mm/sec
        
        % Add FT frame times
        newRow.frameTimes = {ftFrameTimes};
        
        % Add video-related data
        newRow.badVidFrames = {zeros(size(ftFrameTimes))};
        newRow.meanFlow = {meanFlowMags{iTrial}'};
        if exist('meanFlowMagsBallROI', 'var')
            newRow.meanFlowBall = {meanFlowMagsBallROI{iTrial}'};
        end
        rawVidFrameTimes = rawFt(iTrial).rawVidFrameTimes;
        trialVidFrameTimes = rawVidFrameTimes(...
                rawFt(iTrial).startVidFrame:rawFt(iTrial).endVidFrame);
        trialVidFrameTimes = trialVidFrameTimes - rawVidFrameTimes(rawFt(iTrial).startVidFrame - 1);
        newRow.vidFrameTimes = {trialVidFrameTimes};
        
        
        % Append to main table
        ftData = [ftData; newRow];
    end
end

% Save processed data
save(fullfile(outputDir, [expID, '_ficTracData.mat']), 'ftData');

catch ME; rethrow(ME); end

%% Choose threshold for defining movement epochs

moveThresh = 0.06;
flowYLim = [0 0.8];

try
    
% Scan file names in experiment directory to get the expID
imgDataFiles = dir(fullfile(expDir, ['*trial*.tif']));
expID = imgDataFiles(1).name(1:10);

% Load trial and experiment metadata
expMd = readtable(fullfile(outputDir, [expID, '_expMetadata.csv']));
load(fullfile(outputDir, [expID, '_trialMetadata.mat']), 'trialMetadata');

% Load optic flow data
load(fullfile(outputDir, [expID, '_ficTracData.mat']), 'ftData');
meanFlow = ftData.meanFlow;

figure(1);clf;
set(gcf, 'color', [1 1 1])
nTrials = size(ftData, 1);
for iTrial = 1:nTrials
    
    subaxis(nTrials, 1, iTrial, 'S', 0.02, 'mt', 0.02, 'mb', 0.02);

    % Optic flow data to identify flailing
    currFlow = meanFlow{iTrial};
    currFlow(end) = 0;
    flowFrameTimes = ftData.vidFrameTimes{iTrial};
    plotData = repeat_smooth(currFlow, 20, 'dim', 1, 'smwin', 6);
    plotData = plotData - min(plotData);
    plot(flowFrameTimes, plotData, 'color', 'k');
    hold on;
    plot([ftData.frameTimes{iTrial}(1), ftData.frameTimes{iTrial}(end)], [moveThresh, moveThresh],...
            'linewidth', 0.5, 'color', 'r');
    ylim(flowYLim)
    xlim([0 max(trialMetadata.trialDuration)])
    ylabel('Optic flow (flailing)')
    
end

catch ME; rethrow(ME); end

%% Identify flailing events (must run previous section first)
try
    
ftFrameTimes = ftData.frameTimes;
flowFrameTimes = {};
currTrialData = innerjoin(ftData, trialMetadata);
for iTrial = 1:numel(meanFlow)
    
    flowFrameTimes{iTrial} = ftData.vidFrameTimes{iTrial};
    % Warn user if there's a discrepency in the flow frame times based on total trial duration
    discrepVal = abs(flowFrameTimes{iTrial}(end) - currTrialData.trialDuration(iTrial));
    if discrepVal > 0.5
        warning([num2str(discrepVal, 4), ' sec discrepancy in estimated trial duration from flow', ...
            ' frame times for trial #', num2str(iTrial)]);
    end
end

flailingEvents = behaviorEvent('flailing');
flailingEvents = flailingEvents.append_flow_data(expMd.expID{1}, currTrialData.trialNum, ...
        meanFlow, flowFrameTimes, moveThresh);

flailingEvents.export_csv(outputDir, 'fileNamePrefix', expMd.expID{1});

catch ME; rethrow(ME); end

%% Extract ROI data (after creating ROIs with panels_ROI_GUI)

imgDataType = 'reg'; % 'raw' or 'reg'

try
    
% Scan file names in experiment directory to get the expID
imgDataFiles = dir(fullfile(expDir, ['*trial*.tif']));
expID = imgDataFiles(1).name(1:10);

% Find ROI def files
roiDefFiles = dir(fullfile(outputDir, ['roiDefs_trial*.mat']));
roiData = [];
for iFile = 1:numel(roiDefFiles)
    
    % Get current trial number
    trialNum = get_trialNum(roiDefFiles(iFile).name);
    disp(roiDefFiles(iFile).name);
    
    % Load ROI defs
    load(fullfile(outputDir, roiDefFiles(iFile).name), 'roiDefs');
    
    % Load imaging data
    load(fullfile(outputDir, regexprep(roiDefFiles(iFile).name, 'roiDefs', ['imagingData_', ...
            imgDataType])), 'imgData');
        
    % Reshape imaging data into 1D frames
    sz = size(imgData);
    rsImgData = reshape(permute(imgData, [3 4 1 2]), sz(3), sz(4), []); % --> [plane, volume, pixel]
    
    % Create table for current trial's ROI data
    currRoiData = [];
    for iRoi = 1:numel(roiDefs)
        newRow = table({expID}, trialNum, {roiDefs(iRoi).name}, {roiDefs(iRoi).subROIs}, {[]}, ...
                'VariableNames', {'expID', 'trialNum', 'roiName', 'subROIs', 'rawFl'});
        currRoiData = [currRoiData; newRow];
    end
    
    % Loop through and extract data for each subROI
    trialBaselines = [];
    for iRoi = 1:numel(roiDefs)
        currRoi = roiDefs(iRoi);
        rawData = [];
        for iSubRoi = 1:numel(currRoi.subROIs)
            currSubRoi = currRoi.subROIs(iSubRoi);
            currImgData = squeeze(rsImgData(currSubRoi.plane, :, :));	% --> [volume, pixel]
            mask = currSubRoi.mask(:);                                  % --> [pixel]
            currImgData = currImgData(:, mask)';                        % --> [pixel, volume]
            rawData = [rawData; currImgData];                           % --> [pixel, volume]
        end
        
        % Average data across pixels and subROIs
        roiDataAvg = mean(rawData, 1)'; % --> [volume]
        currRoiData.rawFl{iRoi} = roiDataAvg;
        
        % Calculate trial baseline using bottom 5th percentile of whole trial's Fl data
        currDataSorted = sort(roiDataAvg);
        trialBaselines(iRoi, 1) = currDataSorted(round(numel(currDataSorted) * 0.05)); 
        
    end%iRoi   
    
    % Add trial baselines to table
    currRoiData.trialBaseline = trialBaselines;
    
    % Append to main ROI data table
    roiData = [roiData; currRoiData];
    
end%iFile

% Subtract minimum value across entire experiment from all ROI data if it is < 0
expMinVal = min(cell2mat(roiData.rawFl));
if expMinVal < 0
    disp(['Subtracting minimum value of ', num2str(expMinVal) ' from all ROI data'])
    for iRow = 1:size(roiData, 1)
        roiData.rawFl{iRow} = (roiData.rawFl{iRow} - expMinVal) + 1; % Add one to allow division
        roiData.trialBaseline(iRow) = roiData.trialBaseline(iRow) - expMinVal + 1;
    end
end

% Calculate expBaseline for each ROI and add to table
roiList = unique(roiData(:, 'roiName'));
baselineVals = [];
for iRoi = 1:size(roiList, 1)
    currRoiData = innerjoin(roiList(iRoi, :), roiData);
    roiDataSort = sort(cell2mat(currRoiData.rawFl));
    baselineVals(iRoi) = roiDataSort(round(numel(currDataSorted * 0.05)));
end
baselineTable = [roiList, table(baselineVals', 'variableNames', {'expBaseline'})];
roiData = innerjoin(roiData, baselineTable);

% Save processed ROI data
save(fullfile(outputDir, [expID, '_roiData.mat']), 'roiData');

catch ME; rethrow(ME); end


%% Combine two or more ROIs into a new ROI and add it to the ROI def file

combROIs = {'EB', 'BU-L', 'BU-R'};
newROI = 'EB-DAN';

combine_roiDefs(outputDir, combROIs, newROI);

%% Combine PB glomerulus ROIs into EB wedge ROIs
       
glomPairNames = table((1:8)', {'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'}', ...
    {'R1', 'R8', 'R7', 'R6', 'R5', 'R4', 'R3', 'R2'}', 'variablenames', ...
    {'wedge', 'leftGlom', 'rightGlom'});

for iWedge = 1:8
    disp(iWedge)
    combine_roiDefs(outputDir, {glomPairNames.leftGlom{iWedge}, glomPairNames.rightGlom{iWedge}}, ...
            ['EB-', num2str(iWedge)]);
end
disp('Complete')

%% Process odor event data

odorNames = {'EtOH', 'ACV', 'IBA'};
concentrations = {'neat', 'neat', 'neat'};
flowRates = {12, 12, 12};
trialNums = {[1:5], 6, 7};

try 
    
% Load metadata 
expMdFile = dir(fullfile(outputDir, '*expMetadata.csv'));
expID = {expMdFile.name(1:10)};
[expMd, trialMd] = load_metadata(expID, outputDir);

if any(trialMd.usingOptoStim)
    odorStims = odorEvent();
    odorTrials = find(trialMd.usingOptoStim);
    for iCond = 1:numel(odorNames)
        
        % If there are multiple types of stim delivery append them one at a time
        odorName = odorNames{iCond};
        concentration = concentrations{iCond};
        flowRate = flowRates{iCond};
        currTrialNums = trialNums{iCond};
        if isempty(currTrialNums)
            currTrialNums = trialMd.trialNum;
        end
        for iTrial = 1:numel(currTrialNums)
            currTrialMd = trialMd(trialMd.trialNum == currTrialNums(iTrial), :);
            if currTrialMd.usingOptoStim
                stimTiming = currTrialMd.optoStimTiming{:};
                trialDuration = currTrialMd.trialDuration;
                odorStims = odorStims.append_shorthand(expMd.expID{1}, ...
                        currTrialNums(iTrial), stimTiming, trialDuration, odorName, concentration, ...
                        flowRate);
            end
        end%iTrial
    end
    
    % Export to .csv file
    odorStims.export_csv(outputDir, 'fileNamePrefix', expMd.expID{1})
end
catch ME; rethrow(ME); end

%% Copy files over to a grouped analysis data directory
% 
% groupedAnalysisDirName = 'GroupedAnalysisData\new_PPL201_experiments';
groupedAnalysisDirName = 'GroupedAnalysisData_60D05_7f';

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';
analysisDir = fullfile('D:\Dropbox (HMS)\2P Data\Imaging Data', groupedAnalysisDirName);
% expList = {'20201222-1', '20201222-2', '20201228-1', '20201228-2', '20201228-3', '20210102-1', ...
%         '20210102-2', '20210102-3', '20210102-4'};
expList = {'20210118-1'}


for iExp = 1:numel(expList)
    currExpID = expList{iExp};
    disp(currExpID)
    currExpDir = dir(fullfile(parentDir, [currExpID, '*']));
    dataFiles = dir(fullfile(parentDir, currExpDir.name, 'ProcessedData', [currExpID, '*']));
    for iFile = 1:numel(dataFiles)
        copyfile(fullfile(dataFiles(iFile).folder, dataFiles(iFile).name), ...
                fullfile(analysisDir, dataFiles(iFile).name));
    end
end
