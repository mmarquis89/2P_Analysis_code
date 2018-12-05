function normalize_optic_flow(parentDir, sid, varargin)
%=======================================================================================================
% NORMALIZE OPTIC FLOW DATA ACROSS ALL TRIALS IN SESSION
%
% Loads all optic flow data for each trial in a session and normalizes it to the maximum value across
% all trials, then saves normalized flow data for all files in a file named "sid_X_flow_data_norm.mat"
%
% INPUTS:
%       parentDir   = the directory containing the optic flow data files
%
%       sid         = the session ID
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       'OutputDir' = (default: parentDir) the full path to the directory to save the output file in
%========================================================================================================

% Parse optional arguments
p = inputParser;
addParameter(p, 'OutputDir', parentDir);
parse(p, varargin{:});
outputDir = p.Results.OutputDir;

% Create output directory if necessary
if ~isdir(outputDir)
    mkdir(outputDir)
end

% Get list of flow data files
filterStr = ['sid_', num2str(sid), '_tid_*_optic_flow_data.mat'];
flowDataFiles = dir(fullfile(parentDir, filterStr));

% Get list of tids that have flow data
tidList = [];
for iFile = 1:numel(flowDataFiles)
   currFileName = flowDataFiles(iFile).name;
   tidList(iFile) = str2double(regexp(currFileName, '(?<=tid_).*(?=_optic)', 'match'));
end

% Load each file and extract data
allFlowData = []; maxFlowMags = [];
for iFile = 1:numel(flowDataFiles)
    currData = load(fullfile(parentDir, flowDataFiles(iFile).name));
    allFlowData{tidList(iFile)} = currData.meanFlowMag;
    maxFlowMags(tidList(iFile)) = currData.maxFlowMag;
end%iFile

% Normalize flow data 
flyFlowNorm = [];
for iTrial = 1:numel(flowDataFiles)
    flyFlowNorm{tidList(iTrial)} = allFlowData{tidList(iTrial)} ./ max(maxFlowMags);
end%iTrial

% Save file containing normalized flow data for all trials
saveFileName = ['sid_', num2str(sid), '_flow_data_norm.mat'];
save(fullfile(outputDir, saveFileName), 'flyFlowNorm', '-v7.3');

end%function
