function single_plane_heatmap_vid(dataArr, planeNum, infoStruct, cLimRanges, fileName, saveDir, plotTitles, volTitles, figPos, cMapName, frameRate, sigma)
%========================================================================================================================
% CREATE VIDEO OF MULTIPLE dF/F HEATMAP RESPONSES IN A SINGLE PLANE
% 
% Creates a .avi video from a figure broken into a grid of subplots, one for each stimulus type. Each subplot will 
% contain a dF/F heatmap that will be displayed in the video over a number of volumes. The figure can optionally have a 
% title at the top which can change with every volume/frame.
% 
% INPUTS:
%       dataArr     =  a 4-D numeric array with dims [y, x, volume, stimType] containing the dF/F data. The "volume" 
%                      dimension will determine how many frames long the video is.
%
%       planeNum    =  the number of the imaging plane to use for the video
%
%       infoStruct  =  the main imaging data structure containing metadata for the experiment. Specifically, must contain 
%                     the fields "refImg", "MAX_INTENSITY" and "volumeRate".
%
%       cLimRanges  =  an n x 2 vector where each row is in the form [min max] to set the colormap bounds for one plot.
%
%       fileName    =  a string containing the name of the output video file.
%
%       saveDir     =  <OPTIONAL> the directory to save the video file in. Pass [] to generate a popup folder selection 
%                      box instead.
%
%       plotTitles  =  <OPTIONAL> a cell array of strings with length = nStimTypes. Each cell will be the title for a single plot. 
%
%       volTitles  =   <OPTIONAL> a cell array of strings with length == nVols. Each cell will be used as the title for a single 
%                      video frame. 
%
%       figPos      =  <OPTIONAL> position values [x y width height] for the figure window. Pass [] 
%                     to use the default figPos of [50 45 1800 950]
%
%       cMapName    =  <OPTIONAL> name of a colormap function to use in the heatmaps. Pass [] to use 
%                     the default colormap 'bluewhitered'.
%       
%       frameRate   =  <OPTIONAL> the desired frame rate for the output video. Pass [] to use volumeRate as the frame 
%                      rate (resulting in the video playing in real time).
%
%       sigma       =  <OPTIONAL> the degree of gaussian smoothing to be applied to the heatmaps, specified as a 
%                      standard deviation. Pass [] to skip filtering.
%=======================================================================================================================



volumeRate = infoStruct.volumeRate;

% Prompt user for save directory if none was provided
if isempty(saveDir)
    saveDir = uigetdir('B:\Dropbox (HMS)\2P Data\Imaging Data\', 'Select a save directory');
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
if isempty(volTitles)
    volTitles = cell(1, nVols);
end
if isempty(plotTitles)
    plotTitles = cell(1, size(dataArr, 4));
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
    
    % Figure out how many subplots are needed
    nAxes = numSubplots(size(dataArr, 4) + 1);
    
    for iVol = 1:size(dataArr, 3)
        
        % Create figure
        f = figure(652); clf % Choosing arbitrary large number so it doesn't overwrite other figs or the GUI
        f.Position = figPos;
        
        % Plot reference image for the appropriate plane in first subplot
        plotAxes{1} = subaxis(nAxes(1), nAxes(2), 1, 'Spacing', 0.01, 'MB', 0.025);
        imshow(infoStruct.refImg{planeNum}, [0 infoStruct.MAX_INTENSITY]);
        
        % Plot each dF/F heatmap
        for iPlot = 1:size(dataArr, 4)
            
            plotAxes{iPlot + 1} = subaxis(nAxes(1), nAxes(2), iPlot+1, 'Spacing', 0.01, 'MB', 0.025);
            if ~isempty(sigma)
                imagesc(imgaussfilt(dataArr(:, :, iVol, iPlot), sigma));
            else
                imagesc(dataArr(:, :, iVol, iPlot));
            end
            caxis(cLimRanges(iPlot, :));
            colormap(plotAxes{iPlot + 1}, cMapName);
            axis equal; axis off;
            title(plotTitles{iPlot})
            
        end%iPlot
        
        % Add title above all subplots
        suptitle(volTitles{iVol})
        
        % Write frame to video
        writeFrame = getframe(f);
        writeVideo(myVid, writeFrame);
        
    end%iVol
    close(f)
    close(myVid)
end%if
end%function







