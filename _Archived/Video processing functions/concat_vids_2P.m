function concat_vids_2P(vidDir, fileStr, varargin)
%===================================================================================================
% CONCATENATE A SERIES OF VIDEOS INTO A SINGLE FILE
% Concatenates a series of .mp4 videos (e.g. for each trial of a 2P experiment). The new 
% video will be saved in the same folder as the source videos, along with a count of the total
% number of frames in the concatenated video.
% 
% NOTE: the file names will be sorted in alphabetical order to determine sequence of concatenation.
%
% INPUTS:
% 
%       vidDir = the directory containing the videos you want to combine 
%                   e.g. 'U:\2P Behavior Video\2017_07_30\_Movies'
% 
%       fileStr = a filtering string to identify the files to be processed, e.g. '*sid_0_tid*.mp4'. 
%                 Be careful to ensure that only the desired files will meet this specification.
% 
% OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%
%       FrameRate    = (default: 25) the frame rate that the video was acquired at in FPS
%
%       OutputFile   = (default = fileStr + '_AllTrials') the desired name of the output file (minus
%                       the file extension). Default removes all wildcard characters and the fileStr 
%                       and appends '_AllTrials'
%
%       OutputFormat = (default = .mp4) The desired output video file format
%
%===================================================================================================

% Parse optional arguments
p = inputParser;
addParameter(p, 'FrameRate', 25);
addParameter(p, 'OutputFile', [regexprep(regexprep(fileStr, '*', ''), '.mp4', ''), '_AllTrials']);
addParameter(p, 'OutputFormat', '.mp4')
parse(p, varargin{:});
FRAME_RATE = p.Results.FrameRate;
outputFileName = p.Results.OutputFile;
outputFormat = p.Results.OutputFormat;

% Identify files to be processed
vidFiles = dir(fullfile(vidDir, fileStr));
vidNames = sort({vidFiles.name});

% Make sure concatenated video file doesn't already exist
assert(exist(fullfile(vidDir, [outputFileName, outputFormat]), 'file')==0, 'Error: concatenated video file already exists in this directory');

% Create vidWriter
myVidWriter = VideoWriter(fullfile(vidDir, [outputFileName, outputFormat]), 'MPEG-4');
myVidWriter.FrameRate = FRAME_RATE;
open(myVidWriter)
frameCount = 0;

disp('Concatenating videos...')
for iTrial = 1:length(vidNames)
       
    disp(vidNames{iTrial})
            
        % Load movie for the current trial
        myMovie = {};
        myVid = VideoReader(fullfile(vidDir, vidNames{iTrial}));
        while hasFrame(myVid)
            currFrame = readFrame(myVid);
            myMovie(end+1) = {uint8(currFrame)};
        end
        
        % Add frames to movie
        for iFrame = 1:length(myMovie)
            writeVideo(myVidWriter, myMovie{iFrame});
            frameCount = frameCount + 1;
        end
end%iTrial

close(myVidWriter)
clear('myMovie')

% Save frame count log
save(fullfile(vidDir, [outputFileName, '_frameCountLog.mat']), 'frameCount')

disp(['Videos concatenated successfully! Total frames = ', num2str(frameCount)])

end%function