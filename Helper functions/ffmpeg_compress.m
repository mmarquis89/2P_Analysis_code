function ffmpeg_compress(varargin)
%===================================================================================================
% USES FFMPEG TO COMPRESS A VIDEO FILE
% Compresses a target video using a command line call to ffmpeg. Currently only supports H.264 
% encoding (with variable compression quality), but could easily be modified to use other formats
% if necessary.
%
% OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%
%       'SourceVid'  = (default: ask user) the path to the video file to be compressed
%
%       'OutputFileName' = a string to append to the file name of the video. The default value is 
%                          the name of the input video with '_compressed' appended to it
%
%       'CRF' = (default: 23) a number specifying the "constant rate factor" for the compression, to
%               determine the video quality. Smaller numbers are less compression, and values from
%               18-28 are generally reasonable.
%
%===================================================================================================
try
% Parse optional arguments
p = inputParser;
addParameter(p, 'SourceVid', []);
addParameter(p, 'OutputFileName', []);
addParameter(p, 'CRF', 23);
parse(p, varargin{:});
sourceFile = p.Results.SourceVid;
outputFileName = p.Results.OutputFileName;
crf = p.Results.CRF;


% Prompt user for input vid if none was provided
if isempty(sourceFile)
    [fileName, sourceDir] = uigetfile('*');
    sourceFileFormat = regexp(fileName, '\..*', 'match');
    sourceFile = fullfile(sourceDir, fileName);
end

% Create output file name if none was provided
if isempty(outputFileName)
    sourceFileName = regexp(sourceFile, '.*(?=\.)', 'match');
    outputFileName = [sourceFileName{:}, '_compressed', sourceFileFormat{:}];
end

% Create ffmpeg command string
cmdStr = ['ffmpeg -i "', sourceFile, '" -c libx264 -crf ', num2str(crf), ' "', outputFileName, '"']

% Run command
system(cmdStr, '-echo')


end