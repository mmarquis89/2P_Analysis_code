function make_vid(inputDir, sid, tid, varargin)
%===================================================================================================
% CREATE MOVIES FROM .TIF FILES
% Creates a .avi movie from a directory of .tif files captured by the fly behavior camera. Rather 
% than sorting the files alphabetically, they will be sorted under the assumption that they were 
% acquired using FlyCap2 and therefore named with an increasing counter at the end of the file in 
% the following format: 
%    -n.tif 
% Where n is the 4+ digit number representing the cumulative number of frames captured since the 
% camera was plugged in.
%  
% Inputs:
%
%       inputDir = the directory containing the video frames you want to create a movie from.
%                   e.g. 'U:\2P Behavior Video\2017_07_30'
%
%       sid       = the session ID of the video you want to process.
%
%       tid       = the trial ID of the video you want to process.
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       OutputDir        = (default: dirPath) directory to save the output files in
%
%       FrameRate        = (default: 25) the frame rate that the behavior camera was acquiring at.
%
%===================================================================================================
try
% Parse optional arguments
p = inputParser;
addParameter(p, 'OutputDir', inputDir);
addParameter(p, 'FrameRate', 25);
parse(p, varargin{:});
outputDir = p.Results.OutputDir;
frameRate = p.Results.FrameRate;

% Get list of frame files from the trial's directory
frameFiles = dir(fullfile(inputDir, 'fc2*.tif'));
frameFileNames = {frameFiles.name};

% Extract ID numbers from filenames
match_func = @(x) str2double(regexp(x, '(?<=-)\d*(?=[.]tif)', 'match'));
frameFileNums = cellfun(match_func, frameFileNames);

% Fix sorting order to represent actual time of frame acquisition
[~, sortOrder] = sort(frameFileNums);
frameFilesSorted = frameFiles(sortOrder);

% Padding the trial number with leading zeros if necessary to ensure correct filename sorting
trialStr = ['sid_', num2str(sid), '_tid_', pad(num2str(tid), 3, 'left', '0')]; 

% Create save directory if it doesn't already exist
if ~isdir(outputDir)
    mkdir(outputDir);
end

% Make sure this video doesn't already exist
if exist(fullfile(outputDir, [trialStr, '.mp4']), 'file') == 0
    
    % Create video writer object
    outputVid = VideoWriter(fullfile(outputDir, trialStr), 'Motion JPEG AVI');
    outputVid.FrameRate = frameRate;
    open(outputVid)
    
    % Make sure there's at least one image file in this trial's directory
    if ~isempty(frameFilesSorted)
        currFrames = {frameFilesSorted.name}';
        
        % Write each .tif file to video
        for iFrame = 1:length(currFrames)
            
            disp(currFrames{iFrame})
            
            % Read image
            currImg = imread(fullfile(inputDir, currFrames{iFrame}));
            
            % Write frame to video
            writeVideo(outputVid, currImg);
        end
    end%if
    close(outputVid)
else
    disp(['Video already exists...skipping ', trialStr]);
end%if

catch ME
    write_to_log(getReport(ME), mfilename);
end%try
end%function