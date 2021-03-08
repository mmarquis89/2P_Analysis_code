% ==================================================================================================
% 
% This script is for processing imaging data from experiments conducted since the overhaul of my 
% data storage formatting in summer 2020. Each section will run one step of the process, and for the 
% most part they should be run in the order that they occur in the script.
% 
% The key output files that will be used for analysis after all the processing is done are all 
% either .csv files that can be loaded as tables, or .mat files containing tables, and can be 
% easily joined together using unique identifier columns. They are named as follows, replacing 
% <expID> with the unique experiment identifier (e.g. 20210118-1):
% 
%       <expID>_expMetadata.csv         (Unique ID column:  [expID])
%       <expID>_trialMetadata.mat       (Unique ID columns: [expID][trialNum])
%       <expID>_panelsMetadata.mat      (Unique ID columns: [expID][trialNum])
%       <expID>_ficTracData.mat         (Unique ID columns: [expID][trialNum])
%       <expID>_roiData.mat             (Unique ID columns: [expID][trialNum][roiName])
%       
% And sometimes one or more event (e.g. "flailing", "odor", etc.) data files named as follows:
%       <expID>_event_data_<event type> (Unique ID columns: [expID][trialNum][onsetTime])
%       
% Note: all the try...catch blocks in this script exist to allow the sections to be folded up when 
% you've disabled for/if block and section header code folding in Matlab's settings to avoid the 
% annoying automatic unfolding of your code that they often cause. Any errors that occur during 
% execution of the script will be rethrown.
% 
% Updated by MM 1.18.2021
% ==================================================================================================

%% (1) CHOOSE AN EXPERIMENT DIRECTORY
%---------------------------------------------------------------------------------------------------
% Select the directory that all the experimental data is in (the name should contain the Exp ID and 
% name of the experiment). Inside should be all the raw imaging data and Wombat metadata, as well 
% as a subdirectory called "FicTracData".
% --------------------------------------------------------------------------------------------------
startDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';

expDir = uigetdir(startDir, 'Select an experiment directory');
outputDir = fullfile(expDir, 'ProcessedData');

if ~isfolder(outputDir)
    mkdir(outputDir)
end

%% (optional) MAKE ANATOMY STACK
% Skip this step unless you acquired some high-res anatomy stacks at the end of the experiment. If 
% you do have them, the filenames should start with "Stack_" (although you can change this below).
%---------------------------------------------------------------------------------------------------
try
create_anatomy_stack(expDir, 'FileString', 'Stack_*.tif', 'OutputFilePrefix', 'AnatomyStack', ...
        'OutputDir', outputDir);
catch
   disp('Error: anatomy stack creation failed'); 
end

%% (2) PROCESS RAW IMAGING DATA
% Loads raw ScanImage .tif files, extracts the imaging data, discards the flyback frames and 
% re-saves them as .mat files (along with some volume-averaged reference images).
%---------------------------------------------------------------------------------------------------
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

%% (3) PLANE-WISE NoRMCorre MOTION CORRECTION
% Loads the raw imaging data .mat files created previously and runs NoRMCorre motion correction on 
% them, % one imaging plane at a time (this works better than trying to use 3D correction). To 
% improve the performance, it first clips the highest 0.1% of all fluorescence values from each 
% plane and then smooths each frame with a 2D gaussian filter. Afterwards, it saves the data for 
% each trial in a new .mat file and makes more volume-averaged reference images.
%---------------------------------------------------------------------------------------------------
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

%% (4) CREATE METADATA TABLES
% Re-organizes all the assorted metadata and settings info into three tables: 
%       expMetadata (unique identifier column: [expID])
%       trialMetadata (unique identifier columns: [expID][trialNum])
%       panelsMetadata (unique identifier fields: [expID][trialNum])
% The expMetadata table is saved as a .csv file, and the other two are .mat files because they have 
% some columns that contain vectors or structs.
%---------------------------------------------------------------------------------------------------
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
        
        nPanelsFrames = mD.displayRate * mD.trialDuration;
%         if strcmp(mD.panelsMode, 'open loop')
            xPosFunc = mD.xDimPosFun.func;
            yPosFunc = mD.yDimPosFun.func;
            panelsCycleFrames = numel(mD.xDimPosFun.func);
            panelsCycleTime = panelsCycleFrames / double(mD.displayRate);
            xPosRep = repmat(xPosFunc, 1, ceil(mD.trialDuration / panelsCycleTime));
            yPosRep = repmat(yPosFunc, 1, ceil(mD.trialDuration / panelsCycleTime));
            
            panelsPosX = xPosRep(1:nPanelsFrames);
            panelsPosY = yPosRep(1:nPanelsFrames);
%         else
%             xPosFunc = [];
%             yPosFunc = [];
%         end
        
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

% Save metadata files
writetable(expMd, fullfile(outputDir, [expMd.expID{:}, '_expMetadata.csv']));
save(fullfile(outputDir, [expMd.expID{:}, '_trialMetadata.mat']), 'trialMetadata');
save(fullfile(outputDir, [expMd.expID{:}, '_panelsMetadata.mat']), 'panelsMetadata');

catch ME; rethrow(ME); end

%% (5) PROCESS FicTrac DATA
% Extracts the median luminance from each frame of the raw FicTrac videos and uses it to identify 
% the start and end of the trial in both the videos and the FicTrac output data itself (which have 
% slightly different frame rates). Also determines the times in seconds relative to the start of the 
% video of both the FicTrac data "frames" and the video frames. Then, saves the cropped videos to 
% the processed data directory and the raw FicTrac data + frame times in the root experiment 
% directory (to await further processing)
%---------------------------------------------------------------------------------------------------
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

%% (6) DEFINE OPTIC FLOW ROI(s)
% Selects ROIs for optic flow extraction around the fly (to detect flailing or general movement) or 
% the ball (to detect locomotion specifically) in the behavior videos.
% 
% The workflow is currently set up assuming you will do this using the cropped FicTrac videos that 
% were created in the previous section, but if you want to run the whole script start-to-finish 
% without any user input you could technically move this section to be just above "PROCESS RAW 
% IMAGING DATA" and use the raw FicTrac videos for this at the very beginning of the whole process 
% (you'll just have to make sure you save the ROI files in the processed data directory). 
% 
% For a fly ROI, just use the default save name of 'Behavior_Vid_ROI_Data.mat'
% For a ball ROI, save it as 'Behavior_Vid_ROI_Data_Ball.mat' 
%---------------------------------------------------------------------------------------------------
select_video_ROIs(outputDir);

%% (7) CALCULATE OPTIC FLOW IN THE *FLY* ROI 
% Extracts the mean optic flow from the cropped FicTrac videos in the ROI you previously defined 
% around the fly and saves those values in the processed data directory as "flowMags.mat".
% Looks for an ROI file called "Behavior_Vid_ROI_Data.mat" in that output directory.
%---------------------------------------------------------------------------------------------------
try 
    
% Load or define an ROI around the fly
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

%% (optional) CALCULATE OPTIC FLOW IN THE *BALL* ROI
% Only applicable if the fly was mounted on the ball during this experiment. Does the same thing 
% as the previous section, but instead it looks for an ROI file called 
% "Behavior_Vid_ROI_Data_Ball.mat" and saves the results as "flowMags_ballROI.mat" 
%---------------------------------------------------------------------------------------------------
try 
    
% Load or define an ROI around the ball
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

%% (8) ADDITIONAL FicTrac PROCESSING
% Finish processing the FicTrac data by cropping it down to the imaging acquisition period, adding 
% the fly and/or ball optic flow data, and converting units as needed. Compiles the data into a 
% table:
%       ftData (unique identifier columns: [expID][trialNum]) 
% and then saves it in the output directory.
%---------------------------------------------------------------------------------------------------
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
if exist(fullfile(outputDir, 'flowMags.mat'), 'file')
    load(fullfile(outputDir, 'flowMags.mat'), 'meanFlowMags')
end

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
        if exist('meanFlowMags', 'var')
        newRow.meanFlow = {meanFlowMags{iTrial}'};
        else
            newRow.meanFlow = {[]};
        end
        if exist('meanFlowMagsBallROI', 'var')
            newRow.meanFlowBall = {meanFlowMagsBallROI{iTrial}'};
        else
            newRow.meanFlowBall = {[]};
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

%% (9) EXTRACT ROI DATA
% Only run this section after creating ROIs with "panels_ROI_GUI". Extracts imaging data for each 
% ROI and compiles the average fluorescence and baseline values for each ROI into a table:
%       roiData (unique identifier columns: [expID][trialNum][roiName])
% and then saves it as "roiData.mat".
%
% NOTE: for EPG or EB-DAN imaging experiments, you should run one of the next two sections *before* 
% this one to create compound ROIs.

imgDataType = 'reg'; % 'raw' or 'reg' to specify whether to use motion corrected data.
roiDefsFilePrefix = 'roiDefs'; % will look for files with this string followed by "_trial_XXX.mat" 

%---------------------------------------------------------------------------------------------------
try
    
% Scan file names in experiment directory to get the expID
imgDataFiles = dir(fullfile(expDir, ['*trial*.tif']));
expID = imgDataFiles(1).name(1:10);

% Find ROI def files
roiDefFiles = dir(fullfile(outputDir, [roiDefsFilePrefix, '_trial*.mat']));
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

%% (optional) COMBINE TWO OR MORE ROIs
% Combines two or more existing ROIs into a new ROI and appends that to the roiDefs files 
% for each trial in the experiment (see next section for PB-specific ROI combination).

combROIs = {'EB', 'BU-L', 'BU-R'};
newROI = 'EB-DAN';

%---------------------------------------------------------------------------------------------------
combine_roiDefs(outputDir, combROIs, newROI);

%% (optional) COMBINE PB GLOMERULUS ROIs INTO EB WEDGES
% Only applicable for PB imaging experiments. 
% Combines ROIs for matching left and right PB glomeruli into ROIs representing EB wedges and 
% appends them to the roiDefs files for each trial in the experiment.
% 
% It will look specifically for ROIs named "L1", "L2", etc., including the capitalization.
%---------------------------------------------------------------------------------------------------       
try
glomPairNames = table((1:8)', {'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'}', ...
    {'R1', 'R8', 'R7', 'R6', 'R5', 'R4', 'R3', 'R2'}', 'variablenames', ...
    {'wedge', 'leftGlom', 'rightGlom'});

for iWedge = 1:8
    disp(iWedge)
    combine_roiDefs(outputDir, {glomPairNames.leftGlom{iWedge}, glomPairNames.rightGlom{iWedge}}, ...
            ['EB-', num2str(iWedge)]);
end
disp('Complete')
catch ME; rethrow(ME); end

%% (optional) CHOOSE THRESHOLD FOR DEFINING FLAILING EPOCHS
% Only necessary if the fly was *not* mounted on the ball. Loads and plots all the optic 
% flow data for the experiment so you can set a threshold for defining flailing epochs in the next 
% section.

moveThresh = 0.06;
flowYLim = [0 0.8];
 
%---------------------------------------------------------------------------------------------------
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

%% (optional) IDENTIFY FLAILING EVENTS
% Optional: only necessary if the fly was *not* mounted on the ball. Uses the "moveThresh" threshold
% from the previous section to detect and save all flailing events throughout the experiment as a 
% "behaviorEvent" object.
%---------------------------------------------------------------------------------------------------
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

%% (optional) PROCESS ODOR EVENT DATA
% Only applicable for experiments using an odor stim.
% Uses the following information about odor stim identity and delivery to create and save odor event 
% data for the experiment.

odorNames = {'EtOH', 'ACV', 'IBA'};
concentrations = {'neat', 'neat', 'neat'};
flowRates = {12, 12, 12};
trialNums = {[1:5], 6, 7};

%---------------------------------------------------------------------------------------------------
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

%% (optional) COPY FINAL OUTPUT FILES TO AN ANALYSIS DIRECTORY
% Copies all the final processed data files for one or more experiments over to a specified analysis 
% data directory. Adjust the paths and list of expIDs as necessary before running.

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';

% groupedAnalysisDirName = 'GroupedAnalysisData\new_PPL201_experiments';
% groupedAnalysisDirName = 'GroupedAnalysisData_60D05_7f';
groupedAnalysisDirName = 'EB-DAN_GroupedAnalysisData';
analysisDir = fullfile('D:\Dropbox (HMS)\2P Data\Imaging Data', groupedAnalysisDirName);

% expList = {'20210118-2'};

%---------------------------------------------------------------------------------------------------
try
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
catch ME; rethrow(ME); end

