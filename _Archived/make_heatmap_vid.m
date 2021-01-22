function make_heatmap_vid(plotData, infoStruct, cLimRange, fileName, titleStrings, saveDir, figPos, cMapName, frameRate, sigma)
%==============================================================================================================================
% CREATE VIDEO OF HEATMAPPED dF/F RESPONSES IN ALL PLANES
% 
% Creates a .avi video from a figure broken into a grid of subplots, one for each imaging plane. Each subplot will 
% contain a dF/F heatmap that will be displayed in the video over a number of volumes. The figure can optionally have a 
% title at the top which can change with every volume/frame.
% 
% INPUTS:
%       plotData    =  a 4-D numeric array with dims [y, x, plane, volume] containing the dF/F data. The "volume" 
%                      dimension will determine how many frames long the video is.
%
%       infoStruct  =  the main imaging data structure containing metadata for the experiment. Specifically, must contain 
%                      the fields "nPlanes" and "volumeRate".
%
%       cLimRange   =  a two-element vector [min max] to set the colormap bounds.
%
%       fileName    =  a string containing the name of the output video file.
%
%       titleStings =  a cell array of strings with length == nVols. Each cell will be used as the title for a single 
%                      video frame. 
%
%       saveDir     =  <OPTIONAL> the directory to save the video file in. Pass [] to generate a popup folder selection 
%                      box instead.
%
%       figPos      =  <OPTIONAL> position values [x y width height] for the figure window. Pass [] 
%                      to use the default figPos of [50 45 1800 950]
%
%       cMapName    =  <OPTIONAL> name of a colormap function to use in the heatmaps. Pass [] to use 
%                      the default colormap 'bluewhitered'.
%       
%       frameRate   =  <OPTIONAL> the desired frame rate for the output video. Pass [] to use volumeRate as the frame 
%                      rate (resulting in the video playing in real time).
%
%       sigma       =  <OPTIONAL> the degree of gaussian smoothing to be applied to the heatmaps, specified as a 
%                      standard deviation. Pass [] to skip filtering.
%=======================================================================================================================

nPlanes = infoStruct.nPlanes; 
volumeRate = infoStruct.volumeRate;
nVols = size(plotData, 4); 

% Prompt user for save directory if none was provided
if isempty(saveDir)
    saveDir = uigetdir('D:\Dropbox (HMS)\2P Data\Imaging Data\', 'Select a save directory');
    if saveDir == 0
        % Throw error if user canceled without choosing a directory
        error('ERROR: you must select a save directory or provide one as an argument');
    end
else
    % Create save dir if it doesn't already exist
    if ~isdir(saveDir)
        mkdir(saveDir)
    end
end

% Use default values for any other arguments that were not provided
if isempty(figPos)
    figPos = [50 45 1800 950];
end
if isempty(cMapName)
    cMapName = 'bluewhitered';
end
if isempty(frameRate)
    frameRate = round(volumeRate); % Defaults to real time
end
if isempty(titleStrings)
    titleStrings = cell(1, nVols);
end

% Warn user and offer to cancel save if this video will overwrite an existing file
overwrite = 1;
if exist(fullfile(saveDir, [fileName, '.avi']), 'file') ~= 0
    dlgAns = questdlg('Creating this video will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
    if strcmp(dlgAns, 'No')
        overwrite = 0;
        disp('Saving cancelled')
    end
end

if overwrite
    
    % Create video writer
    myVid = VideoWriter(fullfile(saveDir, fileName));
    myVid.FrameRate = frameRate;
    open(myVid)
    
    for iVol = 1:nVols
        
        % Create fig
        f = figure(523); clf % Choosing a weird number so it doesn't overwrite my GUI or other figs
        f.Position = figPos;
        
        % Figure out how many subplots are needed
        nPlots = numSubplots(nPlanes);
        
        for iPlane = nPlanes:-1:1 % Reverse order so planes go from dorsal --> ventral
            
            % Plot dF/F for each plane
            ax = subaxis(nPlots(1), nPlots(2), iPlane, 'Spacing', 0, 'MB', 0.025);
            if ~isempty(sigma)
                imagesc(imgaussfilt(plotData(:, :, iPlane, iVol), sigma));
            else
                imagesc(plotData(:,:, iPlane, iVol));
            end
            caxis(cLimRange);
            colormap(ax, cMapName);
            axis equal; axis off;
            
            % Label postions
            if iPlane == nPlanes
                title('Ventral')
            elseif iPlane == 1
                title('Dorsal')
            end
            
        end%iPlane
        
        % Add title above all subplots
        suptitle(titleStrings{iVol})
        
        % Write frame to video
        writeFrame = getframe(f);
        writeVideo(myVid, writeFrame);
    end%iVol
    close(f);
end%if
close(myVid);
end%function