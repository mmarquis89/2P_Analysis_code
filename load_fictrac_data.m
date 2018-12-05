function outputData = load_fictrac_data(frameCountInfo, varargin)
%===================================================================================================
%
% Original variables are in the units specified in the FicTrac documentation (mostly radians). All
% three of the calculated speed variables are in rad/sec.
%
% INPUTS:
%
%       frameCountInfo: struct containing fields 'goodTrials' and 'nFrames'
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       'Sid' = (default: 0) the session ID of the FicTrac data to be processed
%
%       'ParentDir' = (default: prompt user) the directory containing the FicTrac data files
%
% OUTPUTS:
%
%       outputData = a structure containing the FicTrac data broken into named fields:
%
%               .frameCounter
%               .dRotCam
%               .dRotError
%               .dRotLab
%               .absOrientCam
%               .absOrientLab
%               .intXY
%               .intHD
%               .moveDirLab
%               .moveSpeed
%               .intForwardMove
%               .intSideMove
%               .timeStamp
%               .seqNum
%           and with additional fields:
%               .yawSpeed
%               .fwSpeed
%               .SideSpeed
%               .droppedFrames --> {trial} [frameList]
%               .resets        --> [trial]
%               .badFtTrials   --> [trial list]
%===================================================================================================

% Parse optional arguments
p = inputParser;
addParameter(p, 'Sid', 0);
addParameter(p, 'ParentDir', []);
parse(p, varargin{:});
sid = p.Results.Sid;
parentDir = p.Results.ParentDir;

nFrames = frameCountInfo.nFrames;
goodTrials = frameCountInfo.goodTrials;

% Prompt user to select input directory if none was provided
if isempty(parentDir)
    parentDir = uigetdir('B:\Dropbox (HMS)\2P Data\Behavior Vids', 'Select a directory of FicTrac Data');
end

% Load FicTrac data from each trial into a single array
dataFiles = dir(fullfile(parentDir, '*tid*.dat'));
allData = []; droppedFrames = []; csvData = []; resets = []; badFtTrials = [];
for iTrial = 1:numel(dataFiles)
    
    % Get current trial number from the file name
    tid = str2double(regexp(dataFiles(iTrial).name, '(?<=tid_).*(?=.dat)', 'match'));
    
    % Load data for current trial
    currData = csvread(fullfile(parentDir, dataFiles(iTrial).name)); % --> [frame, var]
    
    % Mark trial as no good if number of FicTrac frames does not match video
    if currData(end, 1) ~= nFrames
        goodTrials(tid) = 0;
        badFtTrials(end + 1) = tid;
    end
    
    if goodTrials(tid)
        % Deal with dropped frames by repeating rows
        currDroppedFrames = [];
        for iFrame = 1:nFrames
            if size(currData, 1) >= iFrame
                if currData(iFrame, 1) ~= iFrame
                    currData = [currData(1:iFrame-1, :); currData(iFrame-1, :); ...
                                    currData(iFrame:end, :)];
                    currData(iFrame, [1 23]) = [iFrame, currData(iFrame, 23) + 1]; % Update frame and sequence nums
                    currDroppedFrames(end + 1) = iFrame;    % --> [trial, framelist]
                end
            else
                currData(iFrame, :) = currData(end, :);
                currDroppedFrames(end + 1) = iFrame;        % --> [trial, framelist]
            end
        end
        
        currDataCSV = [ones(size(currData, 1), 1) * tid, currData];
        
        % Add to full data arrays
        csvData = [csvData; currDataCSV];
        allData(:,:,tid) = currData; % --> [frame, var, trial]
        droppedFrames{tid} = currDroppedFrames;
        resets(tid) = sum(currData(:,23) == 1) - 1;
    end%if
    
end%iTrial

% Save .csv file of concatenated data for all trials
if ~exist(fullfile(parentDir, 'allTrials.csv'), 'file')
    csvwrite(fullfile(parentDir, 'allTrials.csv'), csvData);
end

% Fill any bad trials with nan
allData(:, :, badFtTrials) = nan;

% Split into individual components --> [frame, var, trial] or [frame, trial]
outputData.frameCounter = squeeze(allData(:, 1, :));
outputData.dRotCam = squeeze(allData(:, [2 3 4], :));          % [X, Y, Z]
outputData.dRotError = squeeze(allData(:, 5, :));
outputData.dRotLab = squeeze(allData(:, [6 7 8], :));          % [X, Y, Z]
outputData.absOrientCam = squeeze(allData(:, [9 10 11], :));   % [X, Y, Z]
outputData.absOrientLab = squeeze(allData(:, [12 13 14], :));  % [X, Y, Z]
outputData.intXY = squeeze(allData(:, [15 16], :));            % [X, Y]
outputData.intHD = squeeze(allData(:, 17, :));
outputData.moveDirLab = squeeze(allData(:, 18, :));
outputData.moveSpeed = squeeze(allData(:, 19, :));
outputData.intForwardMove = squeeze(allData(:, 20, :));
outputData.intSideMove = squeeze(allData(:, 21, :));
outputData.timestamp = squeeze(allData(:, 22, :));
outputData.seqNum = squeeze(allData(:, 23, :));

% Create convenient derived variables
outputData.yawSpeed = [zeros(1, size(outputData.intHD, 2)); diff(smoothdata(unwrap(outputData.intHD, [], 1), 1, 'movmean', 5), 1)];
outputData.fwSpeed = [zeros(1, size(outputData.intHD, 2)); diff(smoothdata(outputData.intForwardMove, 1, 'movmean', 5), 1)];
outputData.sideSpeed = [zeros(1, size(outputData.intHD, 2)); diff(smoothdata(outputData.intSideMove, 1, 'movmean', 5), 1)];

outputData.droppedFrames = droppedFrames;
outputData.resets = resets;
outputData.badFtTrials = badFtTrials;

end