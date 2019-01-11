function make_behavior_vids(vidDataDir, vidSaveDir, sid)
try
    
    % Run full behavior vid processing workflow
    addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster
    write_to_log('Starting behavior vid processing...', mfilename)
    
    % Initialize cluster communication
    c = parcluster;
    write_to_log('Cluster communication opened...', mfilename)
    
    
    %     ---------------------------------------------------------------------------------------------
    
    % Identify block vids
    if ~exist(fullfile(vidDataDir, ['sid_', num2str(sid)]), 'dir')
        blockDataDir = vidDataDir;
    else
        blockDataDir = fullfile(vidDataDir, ['sid_', num2str(sid)]);
        
    end
    blockVids = dir(fullfile(blockDataDir, 'fc2_save*.avi'));
    blockVids = blockVids(~contains({blockVids.name}, 'tid')); % In case there are already single trial vids
    nBlocks = numel(blockVids);
    write_to_log(['nBlocks = ', num2str(nBlocks)],mfilename);
    %     if nBlocks > 9
    %             error('ERROR: behavior video will be sorted incorrectly if there are 10 or more blocks')
    %     end
    
    % Remake block vids so that they can be read correctly by Matlab
    memGB = 2;
    timeLimitMin = 30;
    queueName = 'short';
    jobName = 'remake_block_vids';
    remakeJobArr = [];
    for iBlock = 1:nBlocks
        c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
        inputArgs = {vidDataDir, blockVids(iBlock).name, sid, iBlock-1}
        remakeJobArr{end + 1} = c.batch(@remake_block_vid, 0, inputArgs);
    end
    remakeJobArr = wait_for_jobs(remakeJobArr);
    %
    % Identify new block vids
    blockVids = dir(fullfile(blockDataDir, ['block_vid_sid_', num2str(sid), '_bid_*.avi']));
    blockVids = blockVids(~contains({blockVids.name}, 'tid')); % In case there are already single trial vids
    
    % Separate block videos into invididual trials
    write_to_log('Splitting block videos into individual trials...', mfilename)
    memGB = 4;
    timeLimitMin = 60;
    queueName = 'short';
    jobName = 'split_block_vids'
    splitJobArr = [];
    for iBlock = 1:nBlocks
        
        currBid = regexp(blockVids(iBlock).name, '(?<=bid_).*(?=.avi)', 'match');
        write_to_log(['Block ID = ', currBid{:}], mfilename);
        
        % Check whether this block has already been separated into trials
        if isempty(dir(fullfile(vidDataDir, ['*sid_', num2str(sid), '*bid_', currBid{:}, ...
                '*tid*'])))
            
            write_to_log([fullfile(vidDataDir, ['*sid_', num2str(sid), '*bid_', currBid{:}, ...
                '*tid*']), ' does not exist...splitting frames now'])
            
            % Check whether block is closed-loop (and therefore all one trial)
            mdFileStr = ['metadata*sid_', num2str(sid), '_bid_', num2str(str2double(currBid{:})), '.mat'];
            mdFileName = dir(fullfile(vidDataDir, mdFileStr));
            write_to_log(['mdFileName: ', mdFileName.name], mfilename);
            mData = load(fullfile(vidDataDir, mdFileName.name));
            stimType = mData.metaData.stimTypes{1};
            write_to_log(['Stim type: ', stimType], mfilename);
            
            % Start job
            c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
            if regexp(stimType, 'Closed_Loop')
                inputArgs = {vidDataDir, blockVids(iBlock).name, ...
                    'closedLoop', 1}
            else
                inputArgs = {vidDataDir, blockVids(iBlock).name}
            end
            splitJobArr{end + 1} = c.batch(@split_block_vids, 0, inputArgs);
        else
            write_to_log(['Skipping block ', currBid{:}, ...
                ' because single-trial videos ...already exist'], mfilename);
        end%if
    end%for
    splitJobArr = wait_for_jobs(splitJobArr);
    
    %
    %   ------------------------------------------------------------------------------------------------
    
    
    % Update vid files to reflect any changes from imaging data cleanup, then move to output dir
    vid_dir_cleanup(vidDataDir, sid);
    
    
    %-----------------------------------------------------------------------------------------------
    
    % Get list of the raw video data directories and trial IDs
    trialVids = dir(fullfile(vidDataDir, '*sid*tid*.avi'));
    tid = 0; trialVidNames = []; tidList = [];
    for i = 1:numel(trialVids)
        currSid = str2double(regexp(trialVids(i).name, '(?<=sid_).*(?=_bid)', 'match'));
        if currSid == sid
            tid = tid + 1;
            trialVidNames{end+1} = trialVids(i).name;
            tidList(end+1) = tid;
        end
    end
    disp(tidList);
    write_to_log(['numel(tidList) = ', num2str(numel(tidList))], mfilename);
    
    write_to_log('Raw trial vids identified...', mfilename)
    
    %---------------------------------------------------------------------------------------------------
    
    %
    % Rename videos and copy them to the output directory
    for iTrial = 1:numel(trialVidNames)
        disp(iTrial)
        trialStr = ['sid_', num2str(sid), '_tid_', pad(num2str(tidList(iTrial)), 3, 'left', '0'), ...
            '.avi'];
        sourceFileName = fullfile(vidDataDir, trialVidNames{iTrial});
        destFileName = fullfile(vidSaveDir, trialStr);
        copyfile(sourceFileName, destFileName);
    end
    
    %  %     %---------------------------------------------------------------------------------------------------
    
    
    % Count number of frames in the individual trial videos
    allVidFrameCounts = count_vid_frames(vidSaveDir, sid);
    maxFrames = max(allVidFrameCounts);
    
    disp(allVidFrameCounts)
    disp(maxFrames)
    
    write_to_log(['Video frames counted, maxFrames = ', num2str(maxFrames), '...'], mfilename);
    
    
    %---------------------------------------------------------------------------------------------------
    
    %
    % Start job to concatenate processed trial behavior vids
    memGB = 8;
    timeLimitMin = 120;
    queueName = 'short';
    jobName = 'concatRawVids'
    c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
    fileStr = ['*sid_', num2str(sid), '_tid*.avi'];
    outputFileName = ['sid_', num2str(sid), '_AllTrials'];
    inputArgs = {vidSaveDir, fileStr, 'OutputFile', outputFileName}
    concatVidJob = c.batch(@concat_vids, 0, inputArgs);
    
    %
    %     %---------------------------------------------------------------------------------------------------
    %
    %
    %     % % Start jobs to calculate optic flow for all trials
    %     memGB = 1;
    %     b = 0.04;
    %     % if maxFrames <= 2000
    %     %     b = 0.08;
    %     % else
    %     %     b = 0.02;
    %     % end
    %     timeLimitMin = ceil(b * maxFrames);
    %     queueName = 'short';
    %     jobName = 'opticFlowCalc';
    %     c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
    %     roiFilePath = fullfile(vidDataDir, 'Behavior_Vid_ROI_Data.mat');
    %     flowVidDir = fullfile(vidSaveDir, 'opticFlowVids');
    %     flowJobArr = [];
    %     for iTrial = 1:numel(trialVidNames)
    %         disp([vidSaveDir, ' ', num2str(sid), ' ', num2str(tidList(iTrial)), ' ', roiFilePath, ' ', flowVidDir])
    %         inputArgs = {vidSaveDir, sid, tidList(iTrial), roiFilePath, 'OutputDir', flowVidDir};
    %         flowJobArr{iTrial} = c.batch(@single_trial_optic_flow_calc, 0, inputArgs);
    %     end
    %
    %     write_to_log('Optic flow calc jobs submitted...', mfilename)
    %
    %     % Pause execution until all jobs are done
    %     flowJobArr = wait_for_jobs(flowJobArr);
    
    flowCheckComplete = 0;
    while ~flowCheckComplete
        
        c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
        memGB = 1;
        b = 0.04;
        % if maxFrames <= 2000
        %     b = 0.08;
        % else
        %     b = 0.02;
        % end
        timeLimitMin = ceil(b * maxFrames);
        queueName = 'short';
        jobName = 'opticFlowCalc';
        c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
        roiFilePath = fullfile(vidDataDir, 'Behavior_Vid_ROI_Data.mat');
        flowVidDir = fullfile(vidSaveDir, 'opticFlowVids');
        flowJobArr = [];
        
        % Check whether any any of the flow vid data files are missing
        write_to_log('Checking for missing optic flow data...', mfilename);
        missingTids = [];
        for iTrial = 1:numel(trialVidNames)
            if ~exist(fullfile(flowVidDir, ['sid_', num2str(sid), '_tid_', pad(num2str(tidList(iTrial)), 3, 'left', '0'), '_optic_flow_data.mat']), 'file')
                missingTids(end + 1) = tidList(iTrial);
            end
        end
        if ~isempty(missingTids)
            write_to_log(['Missing flow data for the following trials: ', num2str(missingTids)], mfilename);           
            
            for iTrial = 1:numel(missingTids)
                disp([vidSaveDir, ' ', num2str(sid), ' ', num2str(missingTids(iTrial)), ' ', roiFilePath, ' ', flowVidDir])
                inputArgs = {vidSaveDir, sid, missingTids(iTrial), roiFilePath, 'OutputDir', flowVidDir};
                flowJobArr{iTrial} = c.batch(@single_trial_optic_flow_calc, 0, inputArgs);
            end
            
        else
            write_to_log('All optic flow data accounted for', mfilename);
            flowCheckComplete = 1;
        end
        % Pause execution until all jobs are done
        flowJobArr = wait_for_jobs(flowJobArr);
    end
    
    % Normalize optic flow data across all trials
    normalize_optic_flow(flowVidDir, sid, 'OutputDir', vidSaveDir);
    
    write_to_log('Flow data normalized...', mfilename)
    %
    
    % %     %---------------------------------------------------------------------------------------------------
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
    %     for iTrial = 1:numel(trialVidNames)
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
