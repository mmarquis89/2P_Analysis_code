function make_behavior_vids_test(vidDataDir, vidSaveDir, sid)
try
    
    addpath('/home/mjm60/HelperFunctions')
    
    c = parcluster;    
    
    %---------------------------------------------------------------------------------------------  

    
%     imgTypes = {'*.tif', '*.jpg', '*.png', '*.bmp'};
%     vidTypes = {'fcUncompressed.avi', 'fcJPEG.avi', 'fcH264.mp4', 'vdI420.avi', ...
%         'vdUncompressed.avi', 'vdCinepak.avi', 'vdIYUV.avi', 'vdMSVC.avi'};
%     testDir = fullfile(vidDataDir, 'testDir');
%     frame_read_test(testDir, imgTypes, vidTypes)
    
    %---------------------------------------------------------------------------------------------  

    % Move video frames into individual trial directories
    if isdir(fullfile(vidDataDir, 'BlockVids'))
        
        % Identify block directories
        blockDirs = dir(fullfile(vidDataDir, 'BlockVids', ['*sid_', num2str(sid), '*']));
        blockDirs = blockDirs([blockDirs.isdir]);
        nBlocks = numel(blockDirs);
        
        memGB = 4;
        queueName = 'short';
        jobName = 'split_frames'
        splitJobArr = [];
        for iBlock = 1:nBlocks
            
            currBid = regexp(blockDirs(iBlock).name, '(?<=bid_).*', 'match');
            
            % Check whether this block has already been separated into trials
            if isempty(dir(fullfile(vidDataDir, ['*sid_', num2str(sid), '*bid_', currBid{:}, ...
                    '*tid*'])))
                
                % Start job
                nFiles = numel(dir(fullfile(vidDataDir, 'BlockVids', blockDirs(iBlock).name, ...
                    'fc2*.tif')));
                timeLimitMin = ceil(0.02 * nFiles);
                c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
                inputArgs = {vidDataDir, fullfile('BlockVids', blockDirs(iBlock).name)}
                splitJobArr{end + 1} = c.batch(@split_block_frames, 0, inputArgs);
                
            end%if
        end%for
        splitJobArr = wait_for_jobs(splitJobArr);
    end
catch ME
    write_to_log(ME.message, mfilename);
end%try
end%function