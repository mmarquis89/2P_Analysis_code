function extract_sample_frames(dirPath, fileStr, targetPlane)
%===================================================================================================
% EXTRACT SPECIFIC PLANE FROM RAW ANATOMY STACK FILES
%
% Goes through all the anatomy stack files in a directory and extracts the images from a specific 
% plane, saving them as a .tif stack in a new directory. If the file contains more than one channel of imaging 
% data, will save one image from each channel.
%
% INPUTS:
%    dirPath = folder with stack files, e.g. 'B:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_05'
%
%    fileStr = filter string for dir() selection to get anatomy stack files, e.g. '*stack_*'. Be
%              careful to ensure that only the desired files will meet this specification.
%
%    targetPlane = the plane number that you want to extract from each anatomy stack
%===================================================================================================


% Load each .tif stack in a directory, then extract a specific frame from each one and save it as a .png in a new directory
stacks = dir(fullfile(dirPath, fileStr));
nStacks = length(stacks);

% Create directory to save frames in if it does not already exist
savePath = dirPath; %fullfile(dirPath, 'extracted_frames');
if ~isdir(savePath)
    mkdir(savePath)
end

% Extract the target frame from all stacks
for iVol = 1:nStacks
    disp(['Reading ' stacks(iVol).name, '...']);
    currStack = uint16(read_tif(fullfile(dirPath, stacks(iVol).name))); % --> [y, x, plane] OR --> [y, x, plane, channel]
    
    nChannels = size(currStack, 4);
    if nChannels == 1
        % Save sample image
        outputFrame = squeeze(currStack(:,:,targetPlane));
        fileName = ['Plane_', num2str(targetPlane), '_sample_frames.tif'];
        imwrite(outputFrame, fullfile(savePath, fileName), 'WriteMode', 'append', 'Compression', 'none');
    else
        % Save first channel image
        outputFrame = squeeze(currStack(:,:,targetPlane,:)) % --> [y, x, channel]
        fileName = ['Channel_1_Plane_', num2str(targetPlane), '_sample_frames.tif'];
        imwrite(squeeze(outputFrame(:,:,1)), fullfile(savePath, fileName), 'WriteMode', 'append', 'Compression', 'none');
        
        % Save second channel image
        outputFrame = squeeze(currStack(:,:,targetPlane,:)) % --> [y, x, channel]
        fileName = ['Channel_2_Plane_', num2str(targetPlane), '_sample_frames.tif'];
        imwrite(squeeze(outputFrame(:,:,2)), fullfile(savePath, fileName), 'WriteMode', 'append', 'Compression', 'none');
    end
    
end
end