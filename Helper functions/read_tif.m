function [ imgData, tifMetadata ] = read_tif( tifPath )
%===================================================================================================
% Reads a .tif stack of scanimage 2P data into an array.
%
% Uses read_patterned_tifdata() so it doesn't break if scanimage saved a file that is >4GB since
% that isn't allowed in normal .tif files, and the official ScanImage tifReaders are based on .mex
% files so they don't always play nice with the UNIX systems on the O2 cluster.
%
% Input:
%       tifPath = String specifying path to .tif file
%
% Output:
%       imgData     = image data array: [Lines, Pixels, Planes, Volumes]
%
%       tifMetadata = structure containing all the ScanImage metadata for the .tif file
%
% ==================================================================================================

write_to_log(['Attempting to open ', tifPath], mfilename);

% Open .tif file and metadata
if exist(regexprep(tifPath, '.tif', '.mat'), 'file')
    
    % Just load a .mat file for the data, if available
%     write_to_log('Found .mat file to use instead of .tif', mfilename);
    load(regexprep(tifPath, '.tif', '.mat')) % "tifData"
    if exist('tifData', 'var') && ~exist('imgData', 'var')
        imgData = tifData;
    end
    
    % Load just the metadata from the .tif file
    [~, tifMetadata] = read_patterned_tifdata(tifPath, []);
else
    % Otherwise get everything from the .tif file
%     write_to_log('Using .tif file for data and metadata');
    [imgData, tifMetadata] = read_patterned_tifdata(tifPath);
    
    % Save image dimensions
    nLines = size(imgData, 1);
    nPixels = size(imgData, 2);
    
    % Extract info from image description
    try
        frameString = tifMetadata.tifinfo.Software;
    catch
        frameString = tifMetadata.tifinfo.ImageDescription;
    end
    nChannels = length(frameStringKeyLookup(frameString, 'scanimage.SI.hChannels.channelSave'));
    nPlanes = frameStringKeyLookup(frameString, 'hFastZ.numFramesPerVolume');
    
    % Calculate number of volumes per trial
    if nChannels == 1
        nVolumes = size(imgData, 3) / nPlanes;
    else
        nVolumes = size(imgData, 3) / (nPlanes * 2);
    end
    
    % Save image stack to array, accounting for multiple channels if necessary
    if nChannels == 1
        tifDataReshaped = reshape(imgData, [nLines, nPixels, nPlanes, nVolumes]);
        imgData = tifDataReshaped;
    else
        chanData = [];
        for iChannel = 1:nChannels
            chanData(:,:,:,:, iChannel) = reshape(imgData(:,:,iChannel:nChannels:end), [nLines, ...
                nPixels, nPlanes, nVolumes]);
        end
        imgData = chanData;
    end
end

% Close tif file handle
fclose(tifMetadata.fid);


end

