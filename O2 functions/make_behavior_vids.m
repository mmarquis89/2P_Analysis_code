function make_behavior_vids(vidDataDir, vidSaveDir, sid)
try
    
    % Run full behavior vid processing workflow
    addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster
    write_to_log('Starting behavior vid processing...', mfilename)
    
    % Initialize cluster communication
    c = parcluster;
    write_to_log('Cluster communication opened...', mfilename)
    
    
%     ---------------------------------------------------------------------------------------------
   
    
    % Move video frames into individual trial directories
    if isdir(fullfile(vidDataDir, 'BlockVids'))
        
        write_to_log('Splitting video frames into individual trials...', mfilename)
        
        % Identify block directories
        blockDirs = dir(fullfile(vidDataDir, 'BlockVids', ['*sid_', num2str(sid), '*']));
        blockDirs = blockDirs([blockDirs.isdir]);
        nBlocks = numel(blockDirs);
        write_to_log(['nBlocks = ', num2str(nBlocks)],mfilename);
        if nBlocks > 9
            error('ERROR: behavior video will be sorted incorrectly if there are 10 or more blocks')
        end
        
        memGB = 4;
        queueName = 'short';
        jobName = 'split_frames'
        splitJobArr = [];
        for iBlock = 1:nBlocks
            
            if iBlock > 1
               write_to_log('Waiting for previous block job to finish', mfilename); 
               splitJobArr = wait_for_jobs(splitJobArr);
            end
            
            currBid = regexp(blockDirs(iBlock).name, '(?<=bid_).*', 'match');
            write_to_log(['Block ID = ', currBid{:}], mfilename);
            
            % Check whether this block has already been separated into trials
            if isempty(dir(fullfile(vidDataDir, ['*sid_', num2str(sid), '*bid_', currBid{:}, ...
                                '*tid*'])))
                
                write_to_log([fullfile(vidDataDir, ['*sid_', num2str(sid), '*bid_', currBid{:}, ...
                                '*tid*']), ' does not exist...splitting frames now'])
                
                % Check whether block is closed-loop (and therefore all one trial)
                mdFileStr = ['metadata*sid_', num2str(sid), '_bid_', currBid{:}, '.mat'];
                mdFileName = dir(fullfile(vidDataDir, mdFileStr));
                write_to_log(['mdFileName: ', mdFileName.name], mfilename);
                mData = load(fullfile(vidDataDir, mdFileName.name));
                stimType = mData.metaData.stimTypes{1};
                write_to_log(['Stim type: ', stimType], mfilename);
                
                % Start job
                nFiles = numel(dir(fullfile(vidDataDir, 'BlockVids', blockDirs(iBlock).name, ...
                            'fc2*.tif')));
                timeLimitMin = ceil(0.01 * nFiles);
                c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
                if regexp(stimType, 'Closed_Loop')
                    inputArgs = {vidDataDir, fullfile('BlockVids', blockDirs(iBlock).name), ...
                                    'closedLoop', 1}
                else
                    inputArgs = {vidDataDir, fullfile('BlockVids', blockDirs(iBlock).name)}
                end
                splitJobArr{end + 1} = c.batch(@split_block_frames, 0, inputArgs);
            else
                write_to_log(['Skipping block ', currBid{:}, ...
                                ' because single-trial folders ...already exist'], mfilename);
            end%if
        end%for
        splitJobArr = wait_for_jobs(splitJobArr);
    end
    
    
%   ---------------------------------------------------------------------------------------------------
    

    % Update vid files to reflect any changes made during imaging data cleanup
    vid_dir_cleanup(vidDataDir, sid);
    
%     
    %---------------------------------------------------------------------------------------------------

    % Get list of the raw video data directories and trial IDs
    dirContents = dir(fullfile(vidDataDir, '*sid*tid*'));
    tid = 0; dirNames = []; tidList = [];
    for i = 1:numel(dirContents)
        if isdir(fullfile(vidDataDir, dirContents(i).name))
            currSid = str2double(regexp(dirContents(i).name, '(?<=sid_).*(?=_bid)', 'match'));
            if currSid == sid
                tid = tid + 1;
                dirNames{end+1} = fullfile(vidDataDir, [dirContents(i).name]);
                tidList(end+1) = tid;
            end
        end
    end
    disp(tidList);
    write_to_log(['numel(tidList) = ', num2str(numel(tidList))], mfilename);
    
    write_to_log('Raw vid dirs identified...', mfilename)
    
    % Count the number of frames in each directory to find out how long the videos are
    fileCounts = [];
    for iTrial = 1:numel(dirNames)
        currDir = dir(fullfile(dirNames{iTrial}, '*.tif'));
        fileCounts(iTrial) = numel(currDir);
    end
    targetFrames = max(fileCounts);
    
    write_to_log(['Raw vid frames counted (targetFrames = ', num2str(targetFrames), ')...'], mfilename)

    
    %---------------------------------------------------------------------------------------------------
%     
%     
    % Start a job to create each video
    memGB = 1;
    timeLimitMin = ceil(0.035 * targetFrames);
    queueName = 'short';
    jobName = 'makeVid'
    c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
    vidJobArr = [];
    for iTrial = 1:numel(dirNames)
        disp(iTrial);
        inputArgs = {dirNames{iTrial}, sid, tidList(iTrial), 'OutputDir', vidSaveDir}
        vidJobArr{iTrial} = c.batch(@make_vid, 0, inputArgs);
    end
    
    % Pause execution until all jobs are done
    vidJobArr = wait_for_jobs(vidJobArr);
    
    
    %---------------------------------------------------------------------------------------------------
    
    
    % Start job to concatenate behavior vids
    memGB = 8;
    timeLimitMin = 120;
    queueName = 'short';
    jobName = 'concatRawVids'
    c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
    fileStr = ['*sid_', num2str(sid), '_tid*.avi'];
    outputFileName = ['sid_', num2str(sid), '_AllTrials'];
    inputArgs = {vidSaveDir, fileStr, 'OutputFile', outputFileName}
    concatVidJob = c.batch(@concat_vids, 0, inputArgs);
    
    
    %---------------------------------------------------------------------------------------------------
    
    
    % Count number of frames in the individual trial videos
    allVidFrameCounts = count_vid_frames(vidSaveDir, sid);
    maxFrames = max(allVidFrameCounts);
    
    disp(allVidFrameCounts)
    disp(maxFrames)
    
    write_to_log(['Video frames counted, maxFrames = ', num2str(maxFrames), '...'], mfilename);
    
    
    %---------------------------------------------------------------------------------------------------
    
    
    % % Start jobs to calculate optic flow for all trials
    memGB = 1;
%     b = 0.08;
    if maxFrames <= 2000
        b = 0.08;
    else
        b = 0.015;
    end
    timeLimitMin = ceil(b * maxFrames);
    queueName = 'short';
    jobName = 'opticFlowCalc';
    c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
    roiFilePath = fullfile(vidDataDir, 'Behavior_Vid_ROI_Data.mat');
    flowVidDir = fullfile(vidSaveDir, 'opticFlowVids');
    flowJobArr = [];
    for iTrial = 1:numel(dirNames)
        disp([vidSaveDir, ' ', num2str(sid), ' ', num2str(tidList(iTrial)), ' ', roiFilePath, ' ', flowVidDir])
        inputArgs = {vidSaveDir, sid, tidList(iTrial), roiFilePath, 'OutputDir', flowVidDir};
        flowJobArr{iTrial} = c.batch(@single_trial_optic_flow_calc, 0, inputArgs);
        pause(0.5);
    end
    
    % Pause execution until all jobs are done
    flowJobArr = wait_for_jobs(flowJobArr);
    
    write_to_log('Optic flow calc jobs submitted...', mfilename)
    
    % Normalize optic flow data across all trials
    normalize_optic_flow(flowVidDir, sid, 'OutputDir', vidSaveDir);
    
    write_to_log('Flow data calculated and normalized...', mfilename)
%     
%     
%     %---------------------------------------------------------------------------------------------------
%     
%     
%     % Start jobs to create optic flow vids
%     a = 0.006;
%     b = 0.08;
%     memGB = ceil(a * maxFrames);
%     timeLimitMin = ceil(b * maxFrames);
%     queueName = 'short';
%     jobName = 'makeOpticFlowVid'
%     c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
%     frameCountFile = ['sid_', num2str(sid), '_frameCounts.txt'];
%     flowVidDir = fullfile(vidSaveDir, 'opticFlowVids');
%     flowVidJobArr = [];
%     for iTrial = 1:numel(dirNames)
%         inputArgs = {vidSaveDir, sid, tidList(iTrial), frameCountFile, 'OutputDir', flowVidDir};
%         disp([vidSaveDir, ' ', num2str(sid), ' ', num2str(tidList(iTrial)), ' ', frameCountFile, ' ', flowVidDir])
%         flowVidJobArr{iTrial} = c.batch(@create_single_trial_optic_flow_vid, 0, inputArgs);
%     end
%     
%     % Pause execution until all jobs are done
%     flowVidJobArr = wait_for_jobs(flowVidJobArr);
%     
%     write_to_log('Flow vids created...', mfilename)
%     
%     
%     %---------------------------------------------------------------------------------------------------
%     
%     
%     % Start job to concatenate optic flow vids
%     memGB = 1;
%     timeLimitMin = 120;
%     queueName = 'short';
%     jobName = 'concatFlowVids'
%     c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
%     fileStr = ['*sid_', num2str(sid), '_tid*_With_Optic_Flow.avi'];
%     outputFileName = ['sid_', num2str(sid), '_AllTrials_With_Optic_Flow'];
%     inputArgs = {flowVidDir, fileStr, 'OutputFile', outputFileName};
%     concatFlowVidJob = c.batch(@concat_vids, 0, inputArgs);
    
    
    
catch ME
    write_to_log(getReport(ME), mfilename);
end%try
end%function

%===================================================================================================
