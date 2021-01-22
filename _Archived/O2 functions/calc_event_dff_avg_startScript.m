function calc_event_dff_avg_startScript(imgSaveDir, fileStr, sessionDataFile)

% This is just a wrapper for the calc_event_dff_avg() function so that I can dynamically change 
% resource requests based on the data structure size

% Initialize cluster communication
% configCluster;
c = parcluster;

% Extract expDate from directory path
expDate = regexp(imgSaveDir, '(?<=/)20.*(?=/sid)', 'match');
expDate = expDate{:};

% Extract sid from directory path
sid = regexp(sessionDataFile, '(?<=sid_).*(?=\_)', 'match');
sid = sid{:};

% Check array size
m = matfile(fullfile(imgSaveDir, sessionDataFile));
[nColumns, nLines, nPlanes, nVolumes, nTrials] = size(m, 'wholeSession');

% Calculate parameters and start job
memGB = ceil(1.8311e-08 * nColumns * nLines * nPlanes * nVolumes * nTrials);
if memGB > 249
    memGB = 249;
end
timeLimitMin = ceil(3.0518e-08 * nColumns * nLines * nPlanes * nVolumes * nTrials);
if timeLimitMin > 719
    timeLimitMin = 719;
end
queueName = 'short';
jobName = [expDate, '_sid_', num2str(sid), '_calc_event_dff_avg'];
c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
inputArgs = {imgSaveDir, fileStr, sessionDataFile};
c.batch(@calc_event_dff_avg, 0, inputArgs); 
end