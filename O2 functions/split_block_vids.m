function split_block_vids(vidDataDir, blockVidName, varargin)
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
    
    % Strip file extension from block vid name, if present
    blockVidName = regexprep(blockVidName, '.avi', '');
    
    
    write_to_log('Starting extraction', mfilename)
    write_to_log(['Block vid name = ', fullfile(vidDataDir, blockVidName)], mfilename);
    if ~closedLoop
        
        % Extract corner luminance from each frame
        tic
        rawVid = VideoReader(fullfile(vidDataDir, [blockVidName, '.avi']));
        currFrame = []; cornerLum = []; frameSD = [];
        while hasFrame(rawVid)
            currFrame =  readFrame(rawVid);
            currROI = currFrame(end-roiDims(1):end, 1:roiDims(2));
            frameSD(end + 1) = std(double(currFrame(:))); % To watch out for artifact white frames
            cornerLum(end + 1) = mean(currROI(:));
            if ~mod(numel(frameSD), 100)
               disp(num2str(numel(frameSD))); 
            end
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
            write_to_log(ME.message, mfilename);
        end
        write_to_log(['First frames detected'], mfilename)
        
        % Eliminate any frames with major changes to median frame luminance (FlyCap2 artifacts)
        badInds = [];
        for iLoc = 1:numel(firstFrameLocs)
            if (firstFrameLocs(iLoc) <= numel(frameSD)) && (frameSD(firstFrameLocs(iLoc)) < 15)
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
    trialVidName = fullfile([blockVidName, '_tid_', ...
        pad(num2str(trialCount-1), 3, 'left', '0')]);
    trialVid = VideoWriter(fullfile(vidDataDir, trialVidName), 'Motion JPEG AVI');
    trialVid.FrameRate = FRAME_RATE;
    open(trialVid);
    
    while hasFrame(rawVid)
        currFrame = uint8(readFrame(rawVid));
        frameCount = frameCount + 1;
        
        if frameCount == newFrameCounts(trialCount)
            
            % Move onto the next trial
            writeVideo(trialVid, currFrame);
            close(trialVid);
            disp(num2str(trialCount));
            write_to_log(['Wrote video for trial #', num2str(trialCount), ' of ', ...
                num2str(numel(newFrameCounts))], mfilename);
            if hasFrame(rawVid)
                trialCount = trialCount + 1;
                frameCount = 0;
                trialVidName = fullfile([blockVidName, '_tid_', ...
                        pad(num2str(trialCount-1), 3, 'left', '0')]);
                trialVid = VideoWriter(fullfile(vidDataDir, trialVidName), 'Motion JPEG AVI');
                trialVid.FrameRate = FRAME_RATE;
                open(trialVid);
            end
        else
            % Write frame to existing trial
            writeVideo(trialVid, currFrame);
        end
    end
    
catch ME
    write_to_log(getReport(ME), mfilename);
end%function