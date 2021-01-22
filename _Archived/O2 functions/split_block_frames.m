function split_block_frames(parentDir, blockDir, varargin)
%===================================================================================================
% SEPARATE BEHAVIOR VIDEO FRAMES INTO INDIVIDUAL TRIALS
%
% Uses a flash of light in the lower-left corner of the first frame of each trial to divide behavior
% videos into individual trials instead of full blocks. The size of the ROI that is inspected can be
% changed to a different size, specified in pixels. Also does a check to make sure the entire frame
% isn't white since that is an occasional failure mode with Point Grey cameras.
%
% INPUTS:
%       parent          = parent directory for this experiment's behavior video
%
%       blockDir        = the directory (relative to parent) containing all of block's vid frames
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       'roiDims' = (default: [100 50]) the size of the rectangle in the lower-left corner to watch
%
%       'ClosedLoop' = (default: 0)
%
%===================================================================================================

try
    
    addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster
    
    % Initialize cluster communication
    c = parcluster;
    write_to_log('Cluster communication opened...', mfilename)
    
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'roiDims', [100 50]);
    addParameter(p, 'ClosedLoop', 0)
    parse(p, varargin{:});
    roiDims = p.Results.roiDims;
    closedLoop = p.Results.ClosedLoop;
    
    write_to_log('Starting extraction', mfilename)
    
    % Identify all image files in directory
    write_to_log(['Block dir: ', fullfile(parentDir, blockDir)],mfilename);
    frameFiles = dir(fullfile(parentDir, blockDir, 'fc2*.tif'));
    frameFileNames = {frameFiles.name};
    nFrames = length(frameFiles);
    
    write_to_log(['Image files identified...nFrames = ', num2str(nFrames)], mfilename)
    
    % Extract ID numbers from filenames
    match_func = @(x) str2double(regexp(x, '(?<=-)\d*(?=[.]tif)', 'match'));
    frameFileNums = cellfun(match_func, frameFileNames);
    
    % Fix sorting order to represent actual order of frame acquisition
    [~, sortOrder] = sort(frameFileNums);
    frameFilesSorted = frameFiles(sortOrder);
    sortedFrameNames = {frameFilesSorted.name};
    
    write_to_log(['Frame files sorted'], mfilename)
    
    if ~closedLoop
        
        % Start jobs to extract corner luminance from each frame
        write_to_log('Extracting corner luminance values...', mfilename)
        cornerLum = []; frameSD = [];
        tic
        
        % Divide frames up into jobs of <= 1000 frames
        nJobs = ceil(nFrames / 1000);
        allFrameNums = [];
        startFrame = 1;
        while startFrame < nFrames
            endFrame = startFrame + 999;
            if endFrame <= nFrames
                allFrameNums{end + 1} = startFrame:endFrame;
            else
                allFrameNums{end + 1} = startFrame:nFrames;
            end
            startFrame = startFrame + 1000;
        end
        nJobs = numel(allFrameNums);
        
        write_to_log(['Starting ', num2str(nJobs), ' luminance reading jobs...'], mfilename);
        
        memGB = 2;
        queueName = 'short';
        jobName = 'extract_lum';
        timeLimitMin = 30;
        lumJobArr = [];
        for iJob = 1:nJobs
            
            % Start job
            c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
            inputArgs = {fullfile(parentDir, blockDir), sortedFrameNames, allFrameNums{iJob}, ...
                roiDims, iJob};
            lumJobArr{end + 1} = c.batch(@read_frame_luminance, 0, inputArgs);
            pause(1);
        end
        pause(5)
        write_to_log('Pausing for 5 sec after starting luminance jobs...', mfilename);
        wait_for_jobs(lumJobArr);
        
        % Load and concatenate all data
        pause(5);
        frameSD = []; cornerLum = []; frameNums = [];
        lumFiles = dir(fullfile(parentDir, blockDir, '*lumData*'));
        for iFile = 1:numel(lumFiles)
            currLum = load(fullfile(parentDir, blockDir, lumFiles(iFile).name));
            frameSD = [frameSD, currLum.frameSD];
            cornerLum = [cornerLum, currLum.cornerLum];
            frameNums = [frameNums, currLum.frameNums];
        end
        
        write_to_log(['Luminance extracted in ', num2str(toc), ' sec'], mfilename)

% --------------------------------------------------------------------------------------------------        
% for iFile = 1:nFrames
%     if ~mod(iFile, 25)
%         write_to_log(['Reading file ', num2str(iFile)], mfilename);
%         disp(['Reading file ', num2str(iFile)])
%     end
%     
%     
%     tifObj = Tiff(fullfile(parentDir, blockDir, sortedFrameNames{iFile}));
%     currFrame = read(tifOjb
%     if iFile == 1
%         currFrame = imread(fullfile(parentDir, blockDir, sortedFrameNames{iFile}));
%         fSz = size(currFrame);
%         currROI = currFrame(end-roiDims(1):end, 1:roiDims(2));
%         
%     elseif iFile < 100
%         %
%         currFrame = imread(fullfile(parentDir, blockDir, sortedFrameNames{iFile}));
%         currROI = currFrame(end-roiDims(1):end, 1:roiDims(2));
%         frameSD(iFile) = std(double(currFrame(:))); % To watch out for artifact white frames
%         
%     else
%         pixRegion = {[fSz(1) - roiDims(1), fSz(1)], [1, roiDims(2)]};
%         currROI = imread(fullfile(parentDir, blockDir, sortedFrameNames{iFile}), 'PixelRegion', ...
%             pixRegion);
%     end
%     cornerLum(iFile) = mean(currROI(:));
% end
%---------------------------------------------------------------------------------------------------

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
        
        % Eliminate any completely white frames
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
        
        
        % Check whether the first frame of any trials was dropped
        newFrameCounts = frameCounts;
        for iTrial = 1:numel(firstFrameLocs)
            if frameCounts(iTrial) > targetFrames
                % Mark both trials as invalid
                firstFrameLocs = [firstFrameLocs(1:iTrial-1), nan, nan, firstFrameLocs(iTrial+1:end)];
                newFrameCounts = [newFrameCounts(1:iTrial - 1), 0, 0, newFrameCounts(iTrial+1:end)];
            end
        end
        
        write_to_log(['Checked for dropped first frames'], mfilename)
    else
        % If it's a closed loop trial, just move all the frames
        firstFrameLocs = 1;
        newFrameCounts = nFrames;
    end
    disp('Moving frames to new directories...')
    for iTrial = 1:numel(firstFrameLocs)
        
        if ~mod(iTrial, 10)
            disp(num2str(iTrial));
            write_to_log(['Moving frames for trial ', num2str(iTrial)], mfilename)
        end
        
        % Create single trial directory, padding the trial number with leading zeros if necessary to
        % ensure correct filename sorting
        nameSep = strsplit(blockDir, filesep);
        blockName = nameSep{end};
        dirName = fullfile(parentDir, [blockName, '_tid_', pad(num2str(iTrial-1), 3, 'left', '0')]);
        if ~isdir(dirName)
            mkdir(dirName)
        end
        
        % Move video frames for that trial to the new directory
        if ~isnan(firstFrameLocs(iTrial))
            for iFile = firstFrameLocs(iTrial):firstFrameLocs(iTrial) + newFrameCounts(iTrial) - 1
                sourceFile = fullfile(parentDir, blockDir, sortedFrameNames{iFile});
                destDir = [dirName, filesep];
                if ~isdir(destDir)
                    mkdir(destDir)
                end
                cmdStr = ['cp "', sourceFile, '" "', destDir, '"'];
                [status, cmdout] = system(cmdStr);
            end
        end
        
    end%iTrial
catch ME
    write_to_log(getReport(ME), mfilename);
end

end%function









