function O2_preprocessing_test(expDate, sid)

disp(expDate)
addpath('/home/mjm60/HelperFunctions')

% Initialize cluster communication
% configCluster;
c = parcluster; 

% Set directory paths
vidDataDir= ['/n/scratch2/mjm60/', expDate, '/BehaviorVideo'];
vidSaveDir= ['/n/scratch2/mjm60/', expDate, '/sid_', num2str(sid), '/BehaviorVideo'];

% Start behavior vid creation job
memGB = 2;
timeLimitMin = 60;
queueName = 'short';
jobName = 'make_behavior_vids_test'
c = set_job_params(c, queueName, timeLimitMin, memGB, jobName);
inputArgs = {vidDataDir, vidSaveDir, sid}
disp(inputArgs)
% O2 = {c.batch(@make_behavior_vids_test, 0, inputArgs)};
c.batch(@make_behavior_vids_test, 0, inputArgs)

% wait_for_jobs(O2)
% outputData = fetchOutputs(O2{:});
% save(fullfile([vidDataDir, 'outputData.mat']), 'outputData')

end%function