function single_trial_optic_flow_calc(vidSaveDir, sid, tid, roiDataFile, varargin)
%=======================================================================================================
% CALCULATE OPTIC FLOW DATA FOR A SINGLE TRIAL
%
% Uses previously defined behavior vid ROI to calculate optic flow around the fly to make behavioral
% annotation in ANVIL easier. 
%
% INPUTS:
%       vidSaveDir  `   = directory with the individual behavior vids for this session
%
%       sid             = the session ID of the files you are processing
%
%       tid             = the trial ID of the files you are processing
%
%       roiDataFile     = FULL PATH to a file containing an ROI around the fly for the behavior video
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       'OutputDir' = (default: vidSaveDir) the full path to the directory to save the output data in
%
%========================================================================================================
try
% %------------------------------------------
% % Create log file for debugging purposes
% myFile = fopen(fullfile('/home/mjm60/flowlogfiles', ...
%         ['sid_', num2str(sid), '_tid_', pad(num2str(tid), 3, 'left', '0'), ...
%         '_log.txt']), 'a');
% fprintf(myFile, [datestr(datetime), ' Log file opened from within job', '\r\n']);
% %------------------------------------------

addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

% %------------------------------------------
% fprintf(myFile, [datestr(datetime), ' Path added', '\r\n']);
% %------------------------------------------

% Parse optional arguments
p = inputParser;
addParameter(p, 'OutputDir', vidSaveDir);
parse(p, varargin{:});
outputDir = p.Results.OutputDir;

% %------------------------------------------
% fprintf(myFile, [datestr(datetime), ' Arguments parsed', '\r\n']);
% %------------------------------------------

% Create output directory if necessary
if ~isdir(outputDir)
    mkdir(outputDir)
end
disp(outputDir)

% %------------------------------------------
% fprintf(myFile, [datestr(datetime), ' Checked for output dir', '\r\n']);
% %------------------------------------------

% Load ROI data
load(roiDataFile);

% Load input behavior video
vidName = ['sid_', num2str(sid), '_tid_', pad(num2str(tid), 3, 'left', '0'), '.avi'];
trialVid = VideoReader(fullfile(vidSaveDir, vidName));
disp(fullfile(vidSaveDir, vidName))

% %------------------------------------------
% fprintf(myFile, [datestr(datetime), ' Trial vid opened', '\r\n']);
% %------------------------------------------

% Calculate optic flow for each movie frame
opticFlow = opticalFlowFarneback;
frameCount = 0;
meanFlowMag = [];
while hasFrame(trialVid)
    frameCount = frameCount + 1;
    currFrame = readFrame(trialVid);
    currFrame = currFrame(:,:,1);
    
    % Calculate optic flow within each ROI
    currFrameFlowData = estimateFlow(opticFlow, currFrame);
    currFlowFly = currFrameFlowData.Magnitude;
    currFlowFly(~roiData) = nan;
    meanFlowMag(frameCount) = (nanmean(currFlowFly(:)));  
end

% %------------------------------------------
% fprintf(myFile, [datestr(datetime), ' All optic flow calculated', '\r\n']);
% %------------------------------------------

% Calculate max flow value in this trial for later normalization
maxFlowMag = max(meanFlowMag(2:end)); % First frame is artificially high so don't count that

% Save optic flow data
save(fullfile(outputDir, ['sid_', num2str(sid), '_tid_', pad(num2str(tid), 3, 'left', '0'), '_optic_flow_data.mat']), 'meanFlowMag', 'maxFlowMag', '-v7.3');

% %------------------------------------------
% fprintf(myFile, [datestr(datetime), ' Optic flow data saved', '\r\n']);
% fclose(myFile);
% %------------------------------------------

disp('Optic flow calculation complete')
catch ME
    write_to_log(getReport(ME), mfilename);
end%try
end%function




