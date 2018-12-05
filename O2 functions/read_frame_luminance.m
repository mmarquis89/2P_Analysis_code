function read_frame_luminance(parentDir, sortedFrameNames, frameNums, roiDims, jobNum)
%===================================================================================================
% READ CORNER LUMINANCE VALUES FOR A SUBSET OF BEHAVIOR VID FRAMES
%
% Takes a list of video frames, loads them, and saves the overall SD as well as the mean 
% luminance in an ROI in the lower-left corner of the image. 
%
% INPUTS:
%       parentDir
%
%       sortedFrameNames
%
%       frameNums
%
%       roiDims
%
%       jobNum
%
%===================================================================================================

% Read luminances for the desired subset of files
cornerLum = []; frameSD = [];
for iFile = 1:numel(frameNums)
    currFrame = imread(fullfile(parentDir, sortedFrameNames{frameNums(iFile)}));
    currROI = currFrame(end-roiDims(1):end, 1:roiDims(2));
    frameSD(iFile) = std(double(currFrame(:))); % To watch out for artifact white frames
    cornerLum(iFile) = mean(currROI(:));
end
% write_to_log(['Lum job #', num2str(jobNum), ': all luminance extracted. Saving...']);

% Save data
saveFileName = fullfile(parentDir, ['lumData_', pad(num2str(jobNum), 3, 'left', '0'), '.mat']);
save(saveFileName, 'cornerLum', 'frameSD', 'frameNums', 'roiDims', '-v7.3');

write_to_log(['Lum job #', num2str(jobNum), ' complete']);

end