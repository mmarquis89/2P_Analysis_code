function nFrames = count_frames(vidFile)
%===================================================================================================
% Loads a video file and returns the number of frames in it.
% 
% Inputs:
%       vidFile = the path (incl. file name) to the file you want to process
%===================================================================================================

nFrames = 0;
try
    myVid = VideoReader(vidFile);
    while hasFrame(myVid)
        readFrame(myVid);
        nFrames = nFrames + 1;
    end 
catch % Because VideoReader throws an error if you try to read an empty video
end
end