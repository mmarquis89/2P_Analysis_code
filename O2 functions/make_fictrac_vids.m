function make_fictrac_vids(expDate, sid)


addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

% Define directory paths
imgDataDir= ['/n/scratch2/mjm60/', expDate, '/ImagingData'];
vidDataDir= ['/n/scratch2/mjm60/', expDate, '/BehaviorVideo'];
imgSaveDir= ['/n/scratch2/mjm60/', expDate, '/sid_', num2str(sid), '/ImagingData'];
vidSaveDir= ['/n/scratch2/mjm60/', expDate, '/sid_', num2str(sid), '/BehaviorVideo'];

% Initialize cluster communication
configCluster;
c = parcluster; 

% Get list of the raw video data directories and trial IDs
dirContents = dir(fullfile(vidDataDir, '*sid*tid*'));
tid = 0; dirNames = []; tidList = [];
for i = 1:numel(dirContents)
    if isdir(fullfile(vidDataDir, dirContents(i).name))
        currSid = str2double(regexp(dirContents(i).name, '(?<=sid_).*(?=_tid)', 'match'));
        if currSid == sid
            tid = tid + 1;
            dirNames{end+1} = fullfile(vidDataDir, [dirContents(i).name]);
            tidList(end+1) = tid;
        end
    end
end
disp(tidList);

% Count the number of frames in each directory to find out how long the videos are
fileCounts = [];
for iTrial = 1:numel(dirNames)
   currDir = dir(fullfile(dirNames{iTrial}, '*.tif'));
   fileCounts(iTrial) = numel(currDir);
end
targetFrames = max(fileCounts);

% Load FicTrac data
load(fullfile(imgSaveDir, 'analysisMetadata.mat')) % --> 'analysisMetadata' struct
ftData = analysisMetadata.ftData;
frameRate = analysisMetadata.FRAME_RATE;

% Load ROI data
load(fullfile(imgSaveDir, 'ROI_Data_Avg.mat')); % --> 'ROIDataAvg', 'ROIDffAvg' ([volume, trial, ROI])

% Pull out variables from FicTrac data
goodFtData.mmXYdata = ftData.intXY * 4.5;                       % --> [frame, axis, trial] (mm)       
goodFtData.mmSpeedData = ftData.moveSpeed * frameRate * 4.5;   % --> [frame, trial] (mm/sec)
goodFtData.mmFWSpeed = ftData.fwSpeed * frameRate * 4.5;       % --> [frame, trial] (mm/sec)
goodFtData.HD = rad2deg(ftData.intHD);                          % --> [frame, trial] (deg)
goodFtData.dHD = rad2deg(ftData.yawSpeed * frameRate);         % --> [frame, trial] (deg/sec)

goodTrialNums = 1:analysisMetadata.nTrials;
goodTrialNums(~analysisMetadata.goodTrials) = [];

goodFtData.frameCounts = fileCounts;

% Scale fluorescence data so min is zero  
goodFl = ROIDataAvg - min(ROIDataAvg(:)); % --> [frame, trial, ROI]

% Normalize to max during entire experiment for each ROI
goodFtData.goodFlNorm = [];
for iROI = 1:size(ROIDataAvg, 3)
    goodFtData.goodFlNorm(:,:,iROI) = goodFl(:,:,iROI) / max(as_vector(goodFl(:,:,iROI)));
end


% Start a job to create each video
memGB = 2;
timeLimitMin = ceil(0.06 * targetFrames);
queueName = 'short';
jobName = 'makeFicTracVid';
c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
vidJobArr = [];
for iTrial = goodTrialNums
    disp(iTrial);
    inputArgs = {vidSaveDir, imgSaveDir, goodFtData, sid, iTrial, 'OutputDir', vidSaveDir};
    vidJobArr{iTrial} = c.batch(@create_single_trial_fictrac_vid, 0, inputArgs);
end

vidJobArr = wait_for_jobs(vidJobArr);

% Start a job to concatenate the newly created FicTrac vids
memGB = 4;
timeLimitMin = ceil(numel(goodTrialNums) * 1.1);
queueName = 'short';
jobName = 'concatFicTracVids'
c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
fileStr = ['*sid_', num2str(sid), '_tid*_With_FicTrac_Plots.avi'];
outputFileName = ['sid_', num2str(sid), '_AllTrials_With_With_FicTrac_Plots'];
inputArgs = {vidSaveDir, fileStr, 'OutputFile', outputFileName};
concatFlowVidJob = c.batch(@concat_vids, 0, inputArgs);





end%function