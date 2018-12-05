function O2_preprocessing(expDate, sid)

disp(expDate)
whos sid
addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

% cd('../../')
disp(cd)

% Initialize cluster communication
% configCluster;
c = parcluster; 

% Set directory paths
imgDataDir= ['/n/scratch2/mjm60/', expDate, '/ImagingData'];
vidDataDir= ['/n/scratch2/mjm60/', expDate, '/BehaviorVideo'];
imgSaveDir= ['/n/scratch2/mjm60/', expDate, '/sid_', num2str(sid), '/ImagingData'];
vidSaveDir= ['/n/scratch2/mjm60/', expDate, '/sid_', num2str(sid), '/BehaviorVideo'];

% Create save directories if they do not exist
disp(isdir(imgSaveDir))
if ~isdir(imgSaveDir)
    mkdir(imgSaveDir)
end
if ~isdir(vidSaveDir)
    mkdir(vidSaveDir)
end
% 
% % Check for irregularities caused by gapless acquisition
% clean_scanimage_files(imgDataDir, sid);
% 
% % Copy any output files over to the vid directory
% if exist(fullfile(imgDataDir, ['sid_', num2str(sid), '_volume_counts.mat']), 'file')
%     copyfile(fullfile(imgDataDir, ['sid_', num2str(sid), '_volume_counts.mat']), ...
%         fullfile(vidDataDir, ['sid_', num2str(sid), '_volume_counts.mat']));
% end
% if ~isempty(dir(fullfile(imgDataDir, ['sid_', num2str(sid), '_bid_*_single_vol_trials.mat'])))
%     copyfile(fullfile(imgDataDir, ['sid_', num2str(sid), '_bid_*_single_vol_trials.mat']), ...
%         vidDataDir);
% end


% Start behavior vid creation job
memGB = 2;
timeLimitMin = 400;
queueName = 'short';
jobName = 'make_behavior_vids'
c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
inputArgs = {vidDataDir, vidSaveDir, sid}
disp(inputArgs)
c.batch(@make_behavior_vids, 0, inputArgs);


% 
% 
% 
% % Make average fluorescence vids
% memGB = 50;
% timeLimitMin = 60;
% queueName = 'short';
% jobName = 'make_fluorescence_vid'
% c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
% inputArgs = {imgDataDir, sid, 'OutputDir', imgSaveDir}
% disp(inputArgs)
% testJob = c.batch(@create_average_fluorescence_vid, 0, inputArgs);
% % 
% % 
% % 
% % 
% % Create anatomy stack
% memGB = 16;
% timeLimitMin = 30;
% queueName = 'short';
% jobName = 'create_anatomy_stack'
% c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
% inputArgs = {imgDataDir, 'OutputDir', imgSaveDir}
% c.batch(@create_anatomy_stack, 0, inputArgs);
% % 
% % Create a second anatomy stack if the files are present in the imgaging data dir
% myFiles = dir(fullfile(imgDataDir, '*Stack2_*.tif'));
% if ~isempty(myFiles)
%     fileStr = '*Stack2_*';
%     outputFilePrefix = 'Stack2';
%     memGB = 16;
%     timeLimitMin = 30;
%     queueName = 'short';
%     jobName = 'create_anatomy_stack'
%     c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
%     inputArgs = {imgDataDir, 'OutputDir', 'imgSaveDir', 'FileString', fileStr, ...
%                 'OutputFilePrefix', outputFilePrefix}
%     c.batch(@create_anatomy_stack, 0, inputArgs);
% end
% % 
% % 
% % 
% % 
% % Run pre-registration routine
% memGB = 150;
% timeLimitMin = 240;
% queueName = 'short';
% jobName = 'pre_reg_processing'
% c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
% inputArgs = {imgDataDir, sid, expDate, 'OutputDir', imgSaveDir}
% preRegJob{1} = c.batch(@preReg_routine_MM, 0, inputArgs);
% 
% % Pause execution until pre-reg job is finished
% preRegJob = wait_for_jobs(preRegJob);
% % 
% % 
% % 
% % 
% % Figure out how big the wholeSession array is
% m = matfile(regexprep(fullfile(imgSaveDir, ['sid_', num2str(sid), '_sessionFile.mat']), '\', '\\\'));
% [~, ~, nPlanes, nVolumes, nTrials] = size(m, 'wholeSession');
% 
% % Run NorRMCorre registration
% fileName = ['sid_', num2str(sid), '_sessionFile.mat'];
% memGB = ceil(0.0007 * nPlanes * nVolumes * nTrials);
% if memGB > 249
%     memGB = 249;
% end
% timeLimitMin = ceil(0.002 * nPlanes * nVolumes * nTrials);
% queueName = 'short';
% jobName = 'NoRMCorre ';
% c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
% inputArgs = {imgSaveDir, fileName}
% regJob{1} = c.batch(@normcorre_registration, 0, inputArgs);
% 
% % Pause execution until reg job is finished
% regJob = wait_for_jobs(regJob);
% % 
% % % 
% % 
% % % 
% % Calculate and save PCA data 
% memGB = ceil(0.0005 * nTrials * nVolumes * nPlanes);
% timeLimitMin = ceil(0.0025 * nTrials * nVolumes * nPlanes);
% queueName = 'short';
% jobName = 'pca_calc';
% c = set_job_params(c, queueName, timeLimitMin, memGB, jobName); 
% inputArgs = {imgSaveDir, ['rigid_sid_', num2str(sid), '_sessionFile.mat']};
% pcaCalcJob{1} = c.batch(@pca_calc, 0, inputArgs);


end%function