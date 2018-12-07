function split_block_vids(vidDataDir, blockVidName, varargin)
%===================================================================================================
% SEPARATE BEHAVIOR BLOCK VIDEO INTO INDIVIDUAL TRIALS
%
% Uses a flash of light in the lower-left corner of the first frame of each trial to divide behavior
% video for an entire block into individual trials. The size of the ROI that is inspected can be
% changed to a different size, specified in pixels. Also does a check to make sure the entire frame
% isn't white since that is an occasional failure mode with Point Grey cameras.
%
% INPUTS:
%       vidDataDir      = parent directory for this experiment's behavior video
%
%       blockVidName    = the name (minus the .avi extension) of the behavior video to be split
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       'roiDims' = (default: [100 50]) the size of the rectangle in the lower-left corner to watch
%
%       'ClosedLoop' = (default: 0) boolean specifying whether the block is only a single trial
%
%       'FRAME_RATE' = (default: 25) the frame rate that the behavior video was acquired at
%
%===================================================================================================
try
    
    addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster
    
%     % Initialize cluster communication
%     c = parcluster;
%     write_to_log('Cluster communication opened...', mfilename)
    
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'roiDims', [100 50]);
    addParameter(p, 'ClosedLoop', 0);
    addParameter(p, 'FRAME_RATE', 25);
    parse(p, varargin{:});
    roiDims = p.Results.roiDims;
    closedLoop = p.Results.ClosedLoop;
    FRAME_RATE = p.Results.FRAME_RATE;
    
    % Remove file extension if it was included in the block name
    if contains(blockVidName, '.avi')
        blockVidName = blockVidName(1:end - 4);
    end
    
    write_to_log('Starting extraction', mfilename)
    write_to_log(num2str(closedLoop), mfilename);
    if ~closedLoop
        
        % Extract corner luminance from each frame
        write_to_log(fullfile(vidDataDir, [blockVidName, '.avi']), mfilename)
        rawVid = VideoReader(fullfile(vidDataDir, [blockVidName, '.avi']));
        write_to_log('Video reader opened', mfilename);
        currFrame = []; cornerLum = []; frameSD = [];frameCount = 0; frameMed = [];
        tic
        while hasFrame(rawVid)
            frameCount = frameCount + 1;
            if ~mod(frameCount, 100)
               write_to_log(['Reading frame ', num2str(frameCount)], mfilename);
               disp(['Reading frame ', num2str(frameCount)]);
            end
            currFrame =  readFrame(rawVid);
            currROI = currFrame(end-roiDims(1):end, 1:roiDims(2));
            frameSD(end + 1) = std(double(currFrame(:))); % To watch out for artifact white frames
            cornerLum(end + 1) = mean(currROI(:));
            frameMed(end + 1) = median(currFrame(:));
        end
         
        write_to_log(['Luminance extracted in ', num2str(toc), ' sec'], mfilename)
        
        % Detect first frames
        cornerLum = cornerLum - median(cornerLum);
        %         peakThresh = 8 * std(cornerLum);
        peakThresh = 0.95 * max(cornerLum);
        MIN_DIST = 5;
        write_to_log('About to find peaks...', mfilename)
        try
            [~, firstFrameLocs] = findpeaks(cornerLum, 'MinPeakHeight', peakThresh, ...
                'MinPeakDistance', MIN_DIST);
        catch ME
            write_to_log(getReport(ME), mfilename);
        end
        write_to_log(['First frames detected'], mfilename)
        
         % Eliminate any frames with major changes to median frame luminance (FlyCap2 artifacts)
         badInds = [];
         meanFrameMed = mean(frameMed);
         baseSubMed = abs(frameMed - meanFrameMed);
         for iLoc = 1:numel(firstFrameLocs)
             if (firstFrameLocs(iLoc) <= numel(frameMed)) && (baseSubMed(firstFrameLocs(iLoc)) > 50)              
                 badInds = [badInds, iLoc];
             end
         end
         firstFrameLocs(badInds) = [];
         
         % Calculate frame counts for each trial
         firstFrameLocs = [1, firstFrameLocs];
         frameCounts = diff([firstFrameLocs, length(cornerLum) + 1]);
         targetFrames = mode(frameCounts);
         
         write_to_log(['frameCounts calculated'], mfilename)
         
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
         
         write_to_log(['Checked for dropped first frames'], mfilename)
         
    else
        % If it's a closed loop trial, just use the entire video
        firstFrameLocs = 1;
        newFrameCounts = nFrames;
    end
    
    % Re-open video reader for full block video
    rawVid = VideoReader(fullfile(vidDataDir, [blockVidName, '.avi']));
    
    
    % ---------------- Write individual trial videos -----------------------------------------------
    trialCount = 1; frameCount = 0;
    
    % Create video writer for first trial
    trialVidName = [blockVidName, '_tid_', pad(num2str(trialCount-1), 3, 'left', '0')];
    trialVid = VideoWriter(fullfile(vidDataDir, trialVidName), 'Motion JPEG AVI');
    trialVid.FrameRate = FRAME_RATE;
    open(trialVid);
    
    write_to_log(num2str(newFrameCounts), mfilename);
    
    while hasFrame(rawVid)
        currFrame = uint8(readFrame(rawVid));
        frameCount = frameCount + 1;
        
        if frameCount > newFrameCounts(trialCount)
            
            % Move onto the next trial
            close(trialVid);
            write_to_log(['Wrote video for trial #', num2str(trialCount), ' of ', ...
                        num2str(numel(newFrameCounts))], mfilename);
            trialCount = trialCount + 1;
            if trialCount > numel(newFrameCounts)
               break 
            end
            frameCount = 1;
            trialVidName = [blockVidName, '_tid_', pad(num2str(trialCount-1), 3, 'left', '0')];
            trialVid = VideoWriter(fullfile(vidDataDir, trialVidName), 'Motion JPEG AVI');
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
    
catch ME
    write_to_log(getReport(ME), mfilename);
end%function