function frame_read_test(parentDir, imgTypes, vidTypes)

roiDims = [100 50];
% imgTypes = {'*.tif', '*.jpeg', '*.png', '*.bmp'};
% vidTypes = {'fcUncompressed.avi', 'fcJPEG2000.avi', 'fcH264.mp4', ...
%               'vdUncompressed.avi', 'vdCinepak.avi', 'vdIYUV.avi', 'vdMSVC.avi'};

% Test image file reading
for iImg = 1:numel(imgTypes)
    tic
    cornerLum = [];
    currFiles = dir(fullfile(parentDir, imgTypes{iImg}));
    nFiles = numel(currFiles);
    for iFile = 1:nFiles
        currFrame = imread(fullfile(parentDir, currFiles(iFile).name));
        currROI = currFrame(end - roiDims(1):end, 1:roiDims(2));
        cornerLum(iFile) = mean(currROI(:));
    end
    logStr = [num2str(nFiles), ' ', imgTypes{iImg}, ' files read in ', num2str(toc), ' seconds'];
    write_to_log(logStr, mfilename)
    disp(logStr)
end

% Test vid file reading
for iVid = 1:numel(vidTypes)
    try
        tic
        cornerLum = [];
        frameCount = [];
        currFile = dir(fullfile(parentDir, vidTypes{iVid}));
        myVid = VideoReader(fullfile(parentDir, currFile(1).name));
        while hasFrame(myVid)
            currFrame = readFrame(myVid);
            currROI = currFrame(end - roiDims(1):end, 1:roiDims(2));
            cornerLum(iFile) = mean(currROI(:));
            frameCount = frameCount + 1;
        end
        logStr = [num2str(frameCount), ' ', vidTypes{iVid}, ' frames read in ', num2str(toc), ' seconds'];
        write_to_log(logStr, mfilename)
        disp(logStr)
    catch
        errStr = ['Error reading video type ', vidTypes{iVid}];
        write_to_log(errStr, mfilename);
        disp(errStr)
    end
end



