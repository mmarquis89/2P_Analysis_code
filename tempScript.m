expList = [{'20201114-1', '20201117-1', '20201117-2', '20201120-2', '20201201-1', '20201201-2', ...
        '20201201-3', '20201203-1', '20201203-2', '20201210-1'}];%%

expList = [{'20201123-1', '20201123-2', '20201124-1', '20201124-2', '202011124-3'}];%%   
    
 startDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';
 
for iExp = 1:numel(expList)
     
    currExpID = expList{iExp};
    expDir = dir(fullfile(startDir, [currExpID, '_*']));
    expDir = fullfile(startDir, expDir.name);
    outputDir = fullfile(expDir, 'ProcessedData');
    
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

end