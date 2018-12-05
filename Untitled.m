
vidDir = 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2018_12_03_exp_1';
sid = 0;
blockVidName = 'fc2_save_2018-12-03-151812-0000.mp4';
roiDims = [100 50];
FRAME_RATE = 25;
roiData = select_video_ROIs('D:\Dropbox (HMS)\2P Data\Behavior Vids');





blockVids = dir(fullfile(vidDir, '*-0000.mp4'));
allFrameCount = 0;
for iBlock = 1:numel(blockVids)
    
    % Read block behavior video
    disp('Reading block vid...')
    opticFlow = opticalFlowFarneback;
    myVid = VideoReader(fullfile(vidDir, blockVids(iBlock).name));
    firstFrame = readFrame(myVid);
    frameCount = 0; frameSD = []; cornerLum = []; meanFlowMag = [];
    while hasFrame(myVid)
        
        frameCount = frameCount + 1;
        
        % Read next frame
        if ~mod(frameCount, 100)
            disp(['Reading frame ', num2str(frameCount), ' of ', ...
                num2str(FRAME_RATE * myVid.Duration)])
        end
        currFrame = readFrame(myVid);
        currFrame = currFrame(:,:,1);
        
        % Calculate optic flow for each movie frame
        currFrameFlowData = estimateFlow(opticFlow, currFrame);
        currFlowFly = currFrameFlowData.Magnitude;
        currFlowFly(~roiData) = nan;
        meanFlowMag(end + 1) = (nanmean(currFlowFly(:)));
        
        % Extract corner luminance from each frame
        currROI = currFrame(end-roiDims(1):end, 1:roiDims(2));
        frameSD(end + 1) = std(double(currFrame(:))); % To watch out for artifact white frames
        cornerLum(end + 1) = mean(currROI(:));
        
    end%while
    
    allFrameCount = allFrameCount + frameCount;
    
    % Detect first frames
    disp('Detecting first frames...')
    cornerLum = cornerLum - median(cornerLum);
    %         peakThresh = 8 * std(cornerLum);
    peakThresh = 0.95 * max(cornerLum);
    MIN_DIST = 5;
    [~, firstFrameLocs] = findpeaks(cornerLum, 'MinPeakHeight', peakThresh, ...
        'MinPeakDistance', MIN_DIST);
    
    % Eliminate any completely white frames
    badInds = [];
    for iLoc = 1:numel(firstFrameLocs)
        if (firstFrameLocs(iLoc) <= numel(frameSD)) && (frameSD(firstFrameLocs(iLoc)) < 15)
            badInds = [badInds, iLoc];
        end
    end
    firstFrameLocs(badInds) = [];
    
    firstFrameLocs(firstFrameLocs == 2) = [];
    
    % Calculate frame counts for each trial
    firstFrameLocs = [1, firstFrameLocs];
    frameCounts = diff([firstFrameLocs, length(cornerLum) + 1]);
    targetFrames = mode(frameCounts);
    
    % Check whether the first frame of any trials was dropped
    newFrameCounts = frameCounts;
    for iTrial = 1:numel(firstFrameLocs)
        if frameCounts(iTrial) > targetFrames
            % Mark both trials as invalid
            firstFrameLocs = [firstFrameLocs(1:iTrial-1), nan, nan, ...
                firstFrameLocs(iTrial+1:end)];
            newFrameCounts = [newFrameCounts(1:iTrial - 1), 0, 0, ...
                newFrameCounts(iTrial+1:end)];
        end
    end
   
end%iBlock



% Save optic flow data
save(fullfile(vidDir, ['sid_', num2str(sid), '_tid_', pad(num2str(iTrial), 3, 'left', '0'), '_optic_flow_data.mat']), 'meanFlowMag', 'maxFlowMag', '-v7.3');

% Calculate max flow value in this trial for later normalization
maxFlowMag = max(meanFlowMag(2:end)); % First frame is artificially high so don't count that






% ---------------- Write individual trial videos -----------------------------------------------
trialCount = 1; frameCount = 0;

% Create video writer for first trial
trialVidName = fullfile(vidDir, ['sid_', num2str(sid), '_bid', num2str(iBlock), '_tid_', ...
    pad(num2str(trialCount), 3, 'left', '0')]);
trialVid = VideoWriter(fullfile(trialVidName), 'Motion JPEG AVI');
trialVid.FrameRate = FRAME_RATE;
open(trialVid);

splitVidArr = [];
for iFrame = 1:size(vidArr, 3)
    
    frameCount = frameCount + 1;
    
    if frameCount > newFrameCounts(trialCount)
        
        % Move onto the next trial
        disp(['Wrote video for trial #', num2str(trialCount), ' of ', ...
            num2str(numel(newFrameCounts))]);
        trialCount = trialCount + 1;
        frameCount = 1;
        trialVidName = fullfile(vidDir, ['sid_', num2str(sid), '_bid', num2str(iBlock), '_tid_', ...
            pad(num2str(trialCount), 3, 'left', '0')]);
        close(trialVid)
        trialVid = VideoWriter(fullfile(trialVidName), 'Motion JPEG AVI');
        trialVid.FrameRate = FRAME_RATE;
        open(trialVid);
    end
    %         if newFrameCounts(trialCount) > 0
    %             writeVideo(trialVid, currFrame);
    %             splitVidArr{trialCount}(:,:, frameCount) = currFrame;
    %         end
    
    % Write frame to existing trial
    currFrame = uint8(vidArr(:,:,iFrame));
    writeVideo(trialVid, currFrame);
    splitVidArr{trialCount}(:,:, frameCount) = currFrame;
end
close(trialVid)

end%iBlock

frameCount = allFrameCount;
save(fullfile(vidDir, ['sid_', num2str(sid), '_AllTrials_frameCountLog.mat']), 'frameCount', '-v7.3');

%% ---------------- Calculate optic flow -----------------------------------------------

roiData = select_video_ROIs('D:\Dropbox (HMS)\2P Data\Behavior Vids');

for iTrial = 1:numel(splitVidArr)
    
    disp(['Calculating flow for trial ', num2str(iTrial)])
    
    % Calculate optic flow for each movie frame
    opticFlow = opticalFlowFarneback;
    meanFlowMag = [];
    for iFrame = 1:size(splitVidArr{iTrial}, 3)
        currFrameFlowData = estimateFlow(opticFlow, splitVidArr{iTrial}(:,:,iFrame));
        currFlowFly = currFrameFlowData.Magnitude;
        currFlowFly(~roiData) = nan;
        meanFlowMag(iFrame) = (nanmean(currFlowFly(:)));
    end
    
    % Calculate max flow value in this trial for later normalization
    maxFlowMag = max(meanFlowMag(2:end)); % First frame is artificially high so don't count that
    
    % Save optic flow data
    save(fullfile(vidDir, ['sid_', num2str(sid), '_tid_', pad(num2str(iTrial), 3, 'left', '0'), '_optic_flow_data.mat']), 'meanFlowMag', 'maxFlowMag', '-v7.3');
    
    
end

%% Normalize optic flow


% Get list of flow data files
filterStr = ['sid_', num2str(sid), '_tid_*_optic_flow_data.mat'];
flowDataFiles = dir(fullfile(vidDir, filterStr));

% Get list of tids that have flow data
tidList = [];
for iFile = 1:numel(flowDataFiles)
    currFileName = flowDataFiles(iFile).name;
    tidList(iFile) = str2double(regexp(currFileName, '(?<=tid_).*(?=_optic)', 'match'));
end

% Load each file and extract data
allFlowData = []; maxFlowMags = [];
for iFile = 1:numel(flowDataFiles)
    currData = load(fullfile(vidDir, flowDataFiles(iFile).name));
    allFlowData{tidList(iFile)} = currData.meanFlowMag;
    maxFlowMags(tidList(iFile)) = currData.maxFlowMag;
end%iFile

% Normalize flow data
flyFlowNorm = [];
for iTrial = 1:numel(flowDataFiles)
    flyFlowNorm{tidList(iTrial)} = allFlowData{tidList(iTrial)} ./ max(maxFlowMags);
end%iTrial

% Save file containing normalized flow data for all trials
saveFileName = ['sid_', num2str(sid), '_flow_data_norm.mat'];
save(fullfile(vidDir, saveFileName), 'flyFlowNorm', '-v7.3');

%% Get vid frame counts

logFile = fopen(fullfile(vidDir, ['sid_', num2str(sid), '_frameCounts.txt']), 'w');
vidFiles = dir(fullfile(vidDir, '*_tid*.avi'));

frameCounts = [];
for iFile = 1:numel(vidFiles)
    disp(['Counting frames for file ', num2str(iFile)])
    frameCounts(iFile) = count_frames(fullfile(vidDir, vidFiles(iFile).name));
    fprintf(logFile, [num2str(frameCounts(iFile)), ',', vidFiles(iFile).name, '\n']);
end
fclose(logFile);



