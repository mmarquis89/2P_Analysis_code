function O2_preprocessing(expDate, sid)

disp(expDate)
whos sid
addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

% cd('../../')
disp(cd)

% Initialize cluster communication
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


% %---------------------------------------------------------------------------------------------------
% 
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
% 
% 
% %---------------------------------------------------------------------------------------------------
% 
% 
% % Start behavior vid creation job
% memGB = 2;
% timeLimitMin = 400;
% queueName = 'short';
% jobName = [expDate, '_sid_', num2str(sid), '_make_behavior_vids']
% c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
% inputArgs = {vidDataDir, vidSaveDir, sid}
% disp(inputArgs)
% c.batch(@make_behavior_vids, 0, inputArgs);
% 
% 
% % ---------------------------------------------------------------------------------------------------
% 
% 
% % Make average fluorescence vids
% memGB = 50;
% timeLimitMin = 60;
% queueName = 'short';
% jobName = [expDate, '_sid_', num2str(sid), '_make_fluorescence_vid']
% c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
% inputArgs = {imgDataDir, sid, 'OutputDir', imgSaveDir}
% disp(inputArgs)
% testJob = c.batch(@create_average_fluorescence_vid, 0, inputArgs);
% 
% 
% %---------------------------------------------------------------------------------------------------

% 
% % Create anatomy stack
% memGB = 16;
% timeLimitMin = 30;
% queueName = 'short';
% jobName = ['create_anatomy_stack_', expDate, '_sid_', num2str(sid)]
% c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
% inputArgs = {imgDataDir, 'OutputDir', imgSaveDir}
% c.batch(@create_anatomy_stack, 0, inputArgs);
% 
% % Create a second anatomy stack if the files are present in the imgaging data dir
% myFiles = dir(fullfile(imgDataDir, '*Stack2_*.tif'));
% if ~isempty(myFiles)
%     fileStr = '*Stack2_*';
%     outputFilePrefix = 'Stack2';
%     memGB = 16;
%     timeLimitMin = 30;
%     queueName = 'short';
%     jobName = [expDate, '_sid_', num2str(sid), '_create_anatomy_stack'];
%     c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
%     inputArgs = {imgDataDir, 'OutputDir', 'imgSaveDir', 'FileString', fileStr, ...
%                 'OutputFilePrefix', outputFilePrefix}
%     c.batch(@create_anatomy_stack, 0, inputArgs);
% end
% 
% 
% % %---------------------------------------------------------------------------------------------------
% 
% 
% % Run pre-registration routine
% memGB = 50;
% timeLimitMin = 240;
% queueName = 'short';
% jobName = [expDate, '_sid_', num2str(sid), '_pre_reg_processing'];
% c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
% inputArgs = {imgDataDir, sid, expDate, 'OutputDir', imgSaveDir}
% preRegJob{1} = c.batch(@preReg_routine_MM, 0, inputArgs);
% 
% % Pause execution until pre-reg job is finished
% preRegJob = wait_for_jobs(preRegJob);
% 
% 
% %---------------------------------------------------------------------------------------------------
% 
% 
% % Run NorRMCorre registration
% sessionDataFiles = dir(fullfile(imgSaveDir, ['sid_', num2str(sid), '_*_plane*sessionFile.mat']));
% nPlanes = numel(sessionDataFiles);
% for iPlane = 1:nPlanes
%     
%     currFile = sessionDataFiles(iPlane).name;
%     m = matfile(fullfile(imgSaveDir, currFile));
%     [nColumns, nLines, nVolumes] = size(m, 'sessionData');
%     clear m
%     
%     memGB = ceil(0.66e-08 * nColumns * nLines * nVolumes);
%     if memGB > 249
%         memGB = 249;
%     end
%     timeLimitMin = ceil(2.5e-08 * nColumns * nLines * nVolumes);
%     if timeLimitMin > 719
%         timeLimitMin = 719;
%     end
%     queueName = 'short';
%     jobPrefix = regexp(currFile, '.*(?=_sessionFile.mat)', 'match');
%     jobName = [expDate, '_', jobPrefix{:}, '_NoRMCorre'];
%     c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
%     inputArgs = {imgSaveDir, currFile};
%     regJobArr{iPlane} = c.batch(@normcorre_registration, 0, inputArgs);
%     
% end%iPlane
% 
% % Pause execution until reg job is finished
% regJobArr = wait_for_jobs(regJobArr);
% 
% 
% %---------------------------------------------------------------------------------------------------
% 
% % Condense reference images into a single file
% refImgFiles = dir(fullfile(imgSaveDir, ['refImage_reg_sid_', num2str(sid), ...
%         '_chan_2_plane*sessionFile.mat']));
% if isempty(refImgFiles)
%     refImgFiles = dir(fullfile(imgSaveDir, ['refImage_reg_sid_', num2str(sid), ...
%         '_chan_1_plane*sessionFile.mat']));
% end
% nPlanes = numel(refImgFiles);
% write_to_log(['Condensing ', num2str(nPlanes), ' reference image files...'], mfilename);
% refImages = []; regTemplates = []; timePointRefImages = [];
% for iPlane = 1:nPlanes
%     currFile = refImgFiles(iPlane).name;
%     planeNum = regexp(currFile, '(?<=_plane_).*(?=_)', 'match'); planeNum = str2double(planeNum{:});
%     write_to_log(['Loading plane ', num2str(planeNum), '...'], mfilename);
%     load(fullfile(imgSaveDir, currFile)); % --> 'refImg', 'regTemplate', 'timePointRefImg'
%     refImages(:, :, planeNum) = refImg;
%     regTemplates(:, :, planeNum) = regTemplate;
%     timePointRefImages(:, :, :, planeNum) = timePointRefImg; 
% end
% save(fullfile(imgSaveDir, ['sid_', num2str(sid), '_refImages.mat']), 'refImages', 'regTemplates', ...
%         '-v7.3');
% write_to_log('Reference images saved', mfilename);
% % --------------------------------------------------------------------------------------------------


% Calculate and save PCA data
sessionDataFiles = dir(fullfile(imgSaveDir, ['rigid*sid_', num2str(sid), ...
    '_*_plane*sessionFile.mat']));
nPlanes = numel(sessionDataFiles);
write_to_log(['Processing PCA data for ', num2str(nPlanes), ' planes'], mfilename);
for iPlane = 1:nPlanes
    write_to_log(['Plane ', num2str(iPlane)], mfilename);
    currFile = sessionDataFiles(iPlane).name;
    write_to_log(currFile, mfilename);
    m = matfile(fullfile(imgSaveDir, currFile));
    [nColumns, nLines, nVolumes] = size(m, 'sessionData');
    clear m
    framePx = nColumns * nLines;
    totalPx = nColumns * nLines * nVolumes;
    
    write_to_log(whos('framePx'), mfilename);
    
    memGB = ceil(2.75e-08 * totalPx);
    if memGB > 249
        memGB = 249;
    end

    if framePx == 128 * 256
        coeffs = [6.7828e-06 -0.0161 3600];
    elseif framePx == 256 * 512
        coeffs = [1.5532e-06 0.0069 3600];
    elseif framePx == 512 * 1024
        coeffs = [3.2e-05 0.3503 3600];
    else
        write_to_log(['Unknown frame size of ' num2str(framePx), ' pixels'], mfilename);
    end
    timeLimitMin = ceil((coeffs(1) * nVolumes^2 + coeffs(2) * nVolumes + coeffs(3)) / 60)*3;
    if timeLimitMin > 719
        timeLimitMin = 719;
    end    
    write_to_log(['timeLimitMin=', num2str(timeLimitMin), ' memGB=', num2str(memGB)], mfilename);
    queueName = 'short';
    jobPrefix = regexp(currFile, '.*(?=_sessionFile.mat)', 'match');
    jobName = [expDate, '_', jobPrefix{:}, '_PCA_calc'];
    c = set_job_params(c, queueName, timeLimitMin, memGB, jobName); 
    inputArgs = {imgSaveDir, currFile};
    pcaCalcJob{iPlane} = c.batch(@single_plane_pca_calc, 0, inputArgs);
    
end%iPlane

% Pause execution until reg job is finished
pcaCalcJob = wait_for_jobs(pcaCalcJob);

%---------------------------------------------------------------------------------------------------


% Condense and save PCA data into a single file
pcaDataFiles = dir(fullfile(imgSaveDir, ['PCA_rigid_reg_sid_', num2str(sid), ...
        '_chan_2_plane*sessionFile.mat']));
if isempty(pcaDataFiles)
    pcaDataFiles = dir(fullfile(imgSaveDir, ['PCA_rigid_reg_sid_', num2str(sid), ...
        '_chan_1_plane*sessionFile.mat']));
end
nPlanes = numel(pcaDataFiles);
write_to_log(['Condensing PCA data for ', num2str(nPlanes), ' planes'], mfilename)
pcaData = []; explained = [];
for iPlane = 1:nPlanes
    write_to_log(['Plane ', num2str(iPlane)], mfilename);
    currFile = pcaDataFiles(iPlane).name;
    planeNum = regexp(currFile, '(?<=plane_).*(?=_session)', 'match'); planeNum = str2double(planeNum{:});
    load(fullfile(imgSaveDir, currFile));   % --> 'coeffs', 'ex'
    pcaData(:, :, :, planeNum) = coeffs;    % --> [y, x, pc, plane]
    explained(:, iPlane) = ex;              % --> [pc, plane]
end
save(fullfile(imgSaveDir, ['PCA_data_sid_', num2str(sid), '.mat']), 'pcaData', 'explained', ...
        '-v7.3');


end%function