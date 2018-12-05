function calc_event_dff_avg_startScript(imgSaveDir, fileStr, sessionDataFile)

% This is just a wrapper for the calc_event_dff_avg() function so that I can dynamically change 
% resource requests based on the data structure size

% Initialize cluster communication
configCluster;
c = parcluster; 

% Check array size
m = matfile(fullfile(imgSaveDir, sessionDataFile));
[~, ~, nPlanes, nVolumes, nTrials] = size(m, 'wholeSession');

% Calculate parameters and start job
memGB = ceil(0.0006 * nPlanes * nVolumes * nTrials);
timeLimitMin = ceil(0.001 * nPlanes * nVolumes * nTrials);
queueName = 'short';
jobName = 'calc_event_dff_avg';
c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
inputArgs = {imgSaveDir, fileStr, sessionDataFile};
c.batch(@calc_event_dff_avg, 0, inputArgs); 
end