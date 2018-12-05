
vidDir = 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2018_12_03_exp_1';
blockVidName = 'fc2_save_2018-12-03-151812-0000.avi';
roiDims = [100 50];
FRAME_RATE = 25;

% Read block behavior video
myVid = VideoReader(fullfile(vidDir, blockVidName));
firstFrame = readFrame(myVid);
firstFrame(:,:,2:end) = [];
vidArr = firstFrame;
while hasFrame(myVid)
    currFrame = readFrame(myVid);
    vidArr(:,:,end + 1) = currFrame(:,:,1);
end

% Extract corner luminance from each frame
frameSD = []; cornerLum = [];
for iFrame = 1:size(vidArr, 3)
    currROI = vidArr(end-roiDims(1):end, 1:roiDims(2), iFrame);
    frameSD(end + 1) = std(double(as_vector(vidArr(:,:,iFrame)))); % To watch out for artifact white frames
    cornerLum(end + 1) = mean(currROI(:));
end

% Detect first frames
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


% ---------------- Write individual trial videos -----------------------------------------------
trialCount = 1; frameCount = 0;

% Create video writer for first trial
trialVidName = fullfile(vidDir, [blockVidName, '_tid_', ...
    pad(num2str(trialCount-1), 3, 'left', '0')]);
trialVid = VideoWriter(fullfile(trialVidName), 'Motion JPEG AVI');
trialVid.FrameRate = FRAME_RATE;
open(trialVid);

for iFrame = 1:size(vidArr, 3)
    currFrame = uint8(vidArr(:,:,iFrame));
    frameCount = frameCount + 1;
    
    if frameCount >= newFrameCounts(trialCount)
        
        % Move onto the next trial
        disp(['Wrote video for trial #', num2str(trialCount), ' of ', ...
            num2str(numel(newFrameCounts))]);
        trialCount = trialCount + 1;
        frameCount = 0;
        trialVidName = fullfile(vidDir, [blockVidName, '_tid_', ...
            pad(num2str(trialCount-1), 3, 'left', '0')]);
        close(trialVid)
        trialVid = VideoWriter(fullfile(trialVidName), 'Motion JPEG AVI');
        trialVid.FrameRate = FRAME_RATE;
        open(trialVid);
        if newFrameCounts(trialCount) > 0
            writeVideo(trialVid, currFrame);
        end
    else
        % Write frame to existing trial
        writeVideo(trialVid, currFrame);
    end
end

%% ---------------- Calculate optic flow -----------------------------------------------

roiData = select_video_ROIs('D:\Dropbox (HMS)\2P Data\Behavior Vids');

for iFrame = 1:size(vidArr, 3)
    
    
end









