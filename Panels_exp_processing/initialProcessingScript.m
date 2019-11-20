

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\20191119-2_38A11_ChR_60D05_7f';
outputDir = fullfile(parentDir, 'ProcessedData');

if ~isdir(outputDir)
    mkdir(outputDir)
end

%% Make anatomy stack
create_anatomy_stack(parentDir, 'FileString', 'Stack2_*.tif', 'OutputFilePrefix', 'AnatomyStack', ...
        'OutputDir', outputDir);
    
    
%% CONSOLIDATE METADATA AND DAQ DATA FOR ALL TRIALS

% Process metadata files
mdFiles = dir(fullfile(parentDir, '*metadata*trial*.mat'));
expMetadata = [];
for iFile = 1:numel(mdFiles)
       
    % Load the file 
    % Contains struct "mD" with fields [trialSettings], [expID], [expDirName], [expDir], 
    % [trialNum], and [SAMPLING_RATE])
    load(fullfile(parentDir, mdFiles(iFile).name)); 
        
    % Flatten structure so the contents of trialSettings are at the root level
    mD = rmfield(setstructfields(mD, mD.trialSettings), 'trialSettings');

    % Concatenate into master structure
    if iFile == 1
        expMetadata = mD;
    else
        expMetadata(iFile) = mD;
    end

end%iFile

% Process DAQ data files
daqDataFiles = dir(fullfile(parentDir, '*daqData*trial*.mat'));
% expDaqData = struct();
for iFile = 1:numel(daqDataFiles)
       
    % Load the file 
    % Contains variabls "outputData", "trialData", and "ColumnLabels"
    currData = load(fullfile(parentDir, daqDataFiles(iFile).name)); 
    
    % Get the trial number of the current file
    currData.trialNum = str2double(regexp(daqDataFiles(iFile).name, '(?<=trial_)...(?=.mat)', ...
            'match', 'once'));
        
    % Add to master structure
    if iFile == 1
        expDaqData = currData;
    else
        expDaqData(iFile) = currData;
    end
    
end%iFile

% Save in processed data directory
save(fullfile(outputDir, 'metadata.mat'), 'expMetadata', '-v7.3');
save(fullfile(outputDir, 'daqData.mat'), 'expDaqData', '-v7.3');

%% PROCESS FICTRAC DATA

ftDir = fullfile(parentDir, 'FicTracData');

% Identify all FicTrac output files
ftVidFiles = dir(fullfile(ftDir, '*trial*.mp4'));
ftDataFiles = dir(fullfile(ftDir, '*trial*.dat'));
ftFrameLogFiles = dir(fullfile(ftDir, '*vidLogFrames*trial*.txt'));

% Display warning if the number of files is non consistent across types
if numel(unique([numel(ftVidFiles), numel(ftDataFiles), numel(ftFrameLogFiles)])) ~= 1
    errordlg('Warning: inconsistent file counts');
end

% Process each trial
currFtData = [];
for iFile = 1:numel(ftVidFiles)
    
    % Get current trial number
    trialNumStr = regexp(ftVidFiles(iFile).name, '(?<=trial_)...', 'match', 'once');
    trialNum = str2double(trialNumStr);
    disp(['Processing FicTrac data for trial #', num2str(trialNum)]);
    
    % Load output video and extract the median luminance from each frame
    currFile = fullfile(ftDir, ftVidFiles(iFile).name);
    currVid = VideoReader(currFile);
    vidData = []; medLum = [];
    frameCount = 0;
    while hasFrame(currVid)
        frameCount = frameCount + 1;
        if ~mod(frameCount, 500)
            disp(['Reading frame ', num2str(frameCount), '...'])
        end
        currFrame = uint16(readFrame(currVid));
        medLum(frameCount) = median(as_vector(currFrame(:,:,1)));
        vidData(:,:, frameCount) = currFrame(:,:,1);
    end
    
    % Determine which video frames mark the beginning and end of the trial
    lumThresh = 0.9 * max(medLum);
    startVidFrame = find(medLum > lumThresh, 1); % Index of first frame with >90% of max luminance
    endVidFrame = find(medLum > lumThresh, 1, 'last');
    
    % Load vid frame log and identify first and last FicTrac data frames from the trial period
    frameLog = csvread(fullfile(ftDir, ftFrameLogFiles(iFile).name));
    startFtFrame = frameLog(startVidFrame) + 1; % Adding one because original is zero-indexed
    endFtFrame = frameLog(endVidFrame) + 1;
    
    % Load main FicTrac data file and pull out data from within the trial period
    rawFtData = csvread(fullfile(ftDir, ftDataFiles(iFile).name));
    currFtData = rawFtData(startFtFrame:endFtFrame, :);
    currFtData(1, :) = currFtData(1,:) - min(currFtData(1, :)); % Align FrameCount to trial start
    
    % Save processed data files along with luminance values and frame log
    ftData(iFile).trialNum = trialNum;
    ftData(iFile).trialData = currFtData;
    ftData(iFile).frameLog = frameLog;
    ftData(iFile).medLum = medLum;
    
    % Write video data from within the trial period to a new file
    vidData = vidData(:, :, startVidFrame:endVidFrame);
    trialVid = VideoWriter(fullfile(outputDir, ['FicTrac_video_trial_', trialNumStr]), 'MPEG-4');
    open(trialVid)
    for iFrame = 1:size(vidData, 3)
        if ~mod(iFrame, 500)
            disp(['Writing frame ', num2str(iFrame), '...'])
        end
        writeVideo(trialVid, uint8(vidData(:,:, iFrame)));
    end
    close(trialVid)

end%iFile
disp('Video processing complete');

% Save data
save(fullfile(outputDir, 'FicTracData.mat'), 'ftData', '-v7.3');


%% PROCESS IMAGING DATA AND METADATA

% Identify imaging data files
imgDataFiles = dir(fullfile(parentDir, '*trial*.tif'));
for iFile = 1:numel(imgDataFiles)

    disp(['Processing file #', num2str(iFile), ' of ', num2str(numel(imgDataFiles))]);
    
    % Get the trial number of the current file
    trialNum = str2double(regexp(imgDataFiles(iFile).name, '(?<=trial_)...(?=.*\.tif)', ...
            'match', 'once'));
        
    % Load the current file
    [imgData, tifMetadata] = read_tif(fullfile(parentDir, imgDataFiles(iFile).name));
    
    % Parse ScanImage metadata
    siMetadata = parse_scanimage_metadata(tifMetadata);
    siMetadata.trialNum = trialNum;
    if iFile == 1
        imagingMetadata = siMetadata;
    else
        imagingMetadata(iFile) = siMetadata;
    end
    
    % Discard flyback frames from imaging data
    nFlybackFrames = siMetadata.SI.hFastZ.numDiscardFlybackFrames;
    imgData(:, :, (end - (nFlybackFrames - 1)):end, :) = []; % --> [y, x, plane, volume]
    
    % Doing these next two steps because for some reason ScanImage offsets fluorescence data values 
    % by some random amount in different individual .tif files...I've found that this approach fixes 
    % that issue and makes the data directly comparable across trials
    
        % Clip bottom 5% of minimum pixel values per frame
        minFrameVals = sort(as_vector(min(min(imgData))));
        clipThresh = minFrameVals(round(length(minFrameVals) * 0.05));
        imgData(imgData < clipThresh) = clipThresh;

        % Then offset so min value = 1 
        imgData = imgData - min(imgData(:));
        imgData = imgData + 1;    
        
    % Save the raw normalized data in a .mat file
    disp('Saving processed data...')
    save(fullfile(outputDir, ['imagingData_raw_trial_', pad(num2str(trialNum), 3, 'left', '0'), '.mat']), ...
            'imgData', '-v7.3');
            
    % Create and save a file containing reference images for each file
    refImages = mean(imgData, 4); % --> [y, x, plane]
    save(fullfile(outputDir, ['refImages_raw_trial_' pad(num2str(trialNum), 3, 'left', '0'), '.mat']), ...
            'refImages', '-v7.3');
    
        
    % ----- Plane-wise (2D) NoRMCorre motion correction -----
    imgDataReg = zeros(size(imgData));
    regTemplates = [];
    for iPlane = 1:size(imgDataReg, 3)
        
        currData = squeeze(imgData(:,:, iPlane, :)); % --> [y, x, volume]
        
        % Set NoRMCorre options
        options_rigid = NoRMCorreSetParms('d1', size(currData, 1), 'd2', size(currData, 2), ...
            'max_shift', [25, 25], ...
            'init_batch', 100, ...
            'us_fac', 50 ...
            );

        % Run registration
        tic; [planeData, ~, regTemplate, ~] = normcorre_batch(currData, options_rigid); toc
                
        % Add data to full volume arrays
        imgDataReg(:, :, iPlane, :) = planeData;        % --> [y, x, plane, volume]
        regTemplates(:, :, iPlane) = regTemplate;       % --> [y, x, plane]
        
    end%iPlane
        
    % Save registered imaging data
    disp('Saving registered data...');
    imgData = imgDataReg;
    save(fullfile(outputDir, ['imagingData_reg_trial_', pad(num2str(trialNum), 3, 'left', '0'), ...
            '.mat']), 'imgData', '-v7.3')
    
    % Save registered data reference images
    refImages = mean(imgData, 4);    % --> [y, x, plane]
    sz = size(imgData);
    imgData = imgData(:, :, :, 1:(100 * floor(sz(4) / 100)));
    imgDeata = reshape(imgData, sz(1), sz(2), sz(3), 100, []);
    timePointRefImages = squeeze(mean(imgData, 4));  % --> [y, x, plane, section of volumes]
    save(fullfile(outputDir, ['refImages_reg_trial_', pad(num2str(trialNum), 3, 'left', '0'), ...
        '.mat']), 'refImages', 'regTemplates', 'timePointRefImages', '-v7.3');
    
end%iFile

% Save the imaging metadata
save(fullfile(outputDir, 'imagingMetadata.mat'), 'imagingMetadata', '-v7.3');

disp('Processing complete')





