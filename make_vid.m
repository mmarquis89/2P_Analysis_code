function make_vid(dirName, sid, tid, expDate)



saveDir = fullfile('/n/scratch2/mjm60/', expDate, 'Movies');
vidName = ['sid_', sid, '_tid_', tid];


% Create video writer object using MPEG-4 H.264 compression for Anvil compatibility
outputVid = VideoWriter(fullfile(saveDir, [vidName, '.mp4']), 'MPEG-4');
outputVid.FrameRate = 25;
open(outputVid)

% Make sure there's at least one image file in this trial's directory
currFiles = dir(fullfile(dirName, '*.tif'));
if ~isempty(currFiles)
    currFrames = sort({currFiles.name}');
    
    % Write each .tif file to video
    for iFrame = 1:length(currFrames)
        
        % Read image
        currImg = imread(fullfile(dirName, currFrames{iFrame}));
        
        % Write frame to video
        writeVideo(outputVid, currImg);
    end
    
    frameCount = length(currFrames);
else
    
    frameCount = 0;
end%if

% Record number of frames in log
myFile = fopen(fullfile(saveDir, ['sid_', sid, '_frameCountLog.txt']), 'a');
fprintf(myFile, [sid, ',', tid, ',', num2str(frameCount), '\n']);
close(outputVid)
fclose(myFile);
end