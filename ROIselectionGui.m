function ROIselectionGui()
%===================================================================================================
% USE REFERENCE IMAGES TO DRAW AND SAVE ROIs FOR A PARTICULAR EXPERIMENT
% 
% Prompts the user to select a reference images file, and then opens a GUI with all planes plotted
% together in the first tab, as well as each plane plotted larger in its own tab. The user can draw
% any number of ROIs on any of the tabs, then save them when satisfied. To draw an ROI, you must 
% first click on the image you will draw it on ("SELECTED" will appear in the image's title when
% you do this), then click the "Draw ROI" button. To switch to a new ROI, click the "Add new ROI" 
% button. When you are finished click "Save current plots" to save all ROIs.
%
% Can also optionally load and plot PCA data and use that to help select ROIs. Note that if you're 
% selecting multi-plane ROIs in the PCA data the plotted ROIs will disappear if you go back to a 
% previous plane, but the ROI will still be included in the saved ROI data.
%
% The ROI data is saved as a cell array called "ROImetadata" containing a struct for each ROI. Each 
% row of the struct represents one individual sub-ROI that was drawn, containing the following data:
%
%           mask   = a 2D logical array specifying the region of the reference image inside the ROI
%           xi     = the X coordinates of each vertex of the sub-ROI
%           yi     = the Y coordinates of each vertex of the sub-ROI
%           plane  = the number of the imaging plane that the sub-ROI is in
%           color  = the RGB value of the color that the sub-ROI was originally plotted in
%           refImg = the reference image that the sub-ROI was drawn on
%
%===================================================================================================

close all

%++++++++++++++++++++++++ INITIALIZATION TASKS +++++++++++++++++++++++++++++++++++++++++++++++++++++

% Prompt user for ref images data file
[dataFile, pathName, ~] = uigetfile('*.mat', 'Select reference image data file', 'D:\Dropbox (HMS)\2P Data\Imaging Data\');

% Load ref images data file
disp(['Loading ' dataFile, '...'])
imgData = load([pathName, dataFile]); % 1 x nPlanes cell array
% refImages = imgData.refImages;
refImages = imgData.regTemplates;
nPlanes = size(refImages, 3);
disp([dataFile, ' loaded'])


% Create hardcoded parameters
maxInts = [];
for iPlane = 1:nPlanes
    maxInts(iPlane) = max(as_vector(refImages(:, :, iPlane)));
end
maxIntDefault = num2str(max(maxInts));
myData.MAX_INTENSITY = str2double(inputdlg('Enter max intensity value for reference images', '', 1, {maxIntDefault}));
MAX_INTENSITY = myData.MAX_INTENSITY; % To control brightness of ref image plots

% Initialize global variables
myData.ROIs = [];
parentIdxROI = 1;
indexROI = 1;
roiPlaneAxes = [];
selected = 0;
ROIplots = [];
ROItext = [];
pcaTab = [];
pcaData = [];

%----------CONSTRUCTING COMPONENTS----------

% Figure
f = figure('position', [50, 45, 1800, 950], 'Tag', 'figure', 'WindowButtonDownFcn', ...
    {@figure_WindowButtonDownFcn}, 'Name', 'ROI Selection GUI', 'NumberTitle', 'off');

%%%% TOP LEVEL TAB GROUP %%%%
baseTabGroup = uitabgroup(f, 'Units', 'Pixels', 'Position', [5 5 1800 940],...
    'Tag', 'baseTabGroup');

% ROI selection tab
roiTab = uitab(baseTabGroup, 'Title', 'ROI selection', 'Tag', 'roiTab');

% ROI subtab group
roiSubtabGroup = uitabgroup(roiTab, 'Units', 'Normalized', 'Position', [0 0 1 0.99],...
    'CreateFcn', {@roiSubtabGroup_CreateFcn}, 'Tag', 'roiSubtabGroup');

% Change "Units" to normalized so components resize automatically with window
f.Units = 'normalized';
baseTabGroup.Units = 'normalized';
roiSubtabGroup.Units = 'normalized';
pcaTab.Units = 'normalized';

%---------------------------------------------------------------------------------------------------
    function roiSubtabGroup_CreateFcn(src, ~)
        
        % Create tabs for each plane and for combined ref images
        roiSelectTabs = [];
        roiRefImages = [];
        
        % Create subTabs
        for iTab = 1:nPlanes+1
            roiSelectTabs{iTab} = uitab(src);
            if iTab == 1
                
                %-----Plot combined image of all planes-----
                roiSelectTabs{1}.Title = 'All planes';
                roiSelectTabs{1}.Tag = 'All planes';
                
                % Calculate necessary axis dimensions
                subplotDim1 = ceil(sqrt(nPlanes));
                subplotDim2 = floor(sqrt(nPlanes));
                if (subplotDim1 * subplotDim2) < nPlanes
                   subplotDim2 = ceil(sqrt(nPlanes)); 
                end                
                padSize = 0.01;
                rightMargin = 0.12;
                axWidth = (1 - rightMargin - ((subplotDim1 + 1) * padSize)) / subplotDim1;
                axHeight = (1 - ((subplotDim2 + 1) * padSize)) / subplotDim2;
                
                % Place axes in correct locations
                roiRefImgAxes = [];
                for yDim = 1:subplotDim2
                    for xDim = 1:subplotDim1
                        xPos = (xDim * padSize) + ((xDim-1) * axWidth);
                        yPos = (yDim * padSize) + ((yDim-1) * axHeight);
                        roiRefImgAxes{yDim, xDim} = axes(roiSelectTabs{1}, 'Units', 'Normalized', 'Position', ...
                            [xPos, yPos, axWidth, axHeight]);
                    end
                end
                

                % Order handles by z-planes, displayed from L-R, top-bottom
%                 roiRefImgAxes = reshape(flipud(roiRefImgAxes)', 1, nPlanes);
                roiRefImgAxes = reshape(flipud(roiRefImgAxes)', 1, numel(roiRefImgAxes));
                if numel(roiRefImgAxes) > nPlanes
                    for iAx = (nPlanes + 1):numel(roiRefImgAxes)
                        axes(roiRefImgAxes{iAx});
                        axis off
                    end
                    roiRefImgAxes = roiRefImgAxes(1:nPlanes);
                end
                
                % Plot reference frames for each plane
                myImages = [];
                for iPlane = 1:nPlanes
                    currAxes = roiRefImgAxes{iPlane};
                    axes(currAxes);
                    axis image; hold on
                    myImages{iPlane} = imshow(refImages(:, :, iPlane), [0 MAX_INTENSITY], ...
                        'InitialMagnification', 'fit', 'Parent', currAxes);
                    myImages{iPlane}.ButtonDownFcn = {@image_ButtonDownFcn};
                    currAxes.Title.String = ['Plane #' num2str(iPlane)];
                    currAxes.UserData = iPlane;
                    axis image
                end
            else
                %-----Plot larger individual plane reference images-----
                roiSelectTabs{iTab}.Title = ['Plane #', num2str(iTab-1)];
                roiPlaneAxes{iTab-1} = axes(roiSelectTabs{iTab}, 'Position', [0.05 0.1 0.8 0.8]);
                axis image; hold on
                roiRefImages{iTab-1} = imshow(refImages(:, :, iTab-1), [0 MAX_INTENSITY], ...
                    'InitialMagnification', 'fit', 'Parent', roiPlaneAxes{iTab-1});
                roiRefImages{iTab-1}.ButtonDownFcn = {@image_ButtonDownFcn};
                roiPlaneAxes{iTab-1}.Title.String = ['Plane #' num2str(iTab-1)];
                roiPlaneAxes{iTab-1}.Tag = 'roiAxes';
                roiPlaneAxes{iTab-1}.UserData = iTab-1;
            end
        end
        
        % Create PCA screening subtab
        pcaTab = uitab(src, 'Title', 'PCA screening', 'Tag', 'pcaTab', 'CreateFcn', ...
                {@pcaTab_CreateFcn});
        
        % Create DrawROI button
        drawROIButton = uicontrol(roiTab, 'Style', 'pushbutton', 'String', 'Draw ROI', 'Units', 'Normalized', ...
            'Position', [0.89, 0.6, 0.1, 0.05], 'FontSize', 12, 'Callback', {@drawROIButton_Callback}, ...
            'Tag', 'drawROIButton');
        
        % Create addROI button
        addROIButton = uicontrol(roiTab, 'Style', 'pushbutton', 'String', 'Add new ROI', 'Units', 'Normalized', ...
            'Position', [0.89, 0.5, 0.1, 0.05], 'FontSize', 12, 'Callback', {@addROIButton_Callback}, ...
            'Tag', 'addROIButton');
        
        % Create clearROIs button
        clearROIButton = uicontrol(roiTab, 'Style', 'pushbutton', 'String', 'Clear ROIs', 'Units', 'Normalized', ...
            'Position', [0.89, 0.4, 0.1, 0.05], 'FontSize', 12, 'Callback', {@clearROIButton_Callback}, ...
            'Tag', 'clearROIButton');
        
        % Create saveROIs button
        ROISaveButton = uicontrol(roiTab, 'Style', 'pushbutton', 'String', 'Save current ROIs', ...
            'Units', 'Normalized', 'Position', [0.89, 0.3, 0.1, 0.05], 'FontSize', 12, ...
            'Callback', {@ROISaveButton_Callback}, 'Tag', 'ROISaveButton');
        
    end

%---------------------------------------------------------------------------------------------------
    function pcaTab_CreateFcn(src, ~)
        
        % Create loadPCA button
        loadPCAButton = uicontrol(roiTab, 'Style', 'pushbutton', ...
                                            'String', 'Load PCA data', ...
                                            'Units', 'Normalized', ...
                                            'Position', [0.89, 0.9, 0.1, 0.05], ...
                                            'FontSize', 12, ...
                                            'Callback', {@loadPCAButton_Callback}, ...
                                            'Tag', 'loadPCAButton');
        % Create plane number slider
        planeNumSlider = uicontrol(roiTab, 'Style', 'slider', ...
                                            'Units', 'normalized', ...
                                            'Position', [0.89, 0.85, 0.1, 0.025], ...
                                            'Callback', {@planeNumSlider_Callback}, ...
                                            'Tag', 'planeNumSlider', ...
                                            'Min', 1, ...
                                            'Max', nPlanes, ...
                                            'SliderStep', [1/nPlanes, 1/nPlanes], ...
                                            'Enable', 'off');
        
    end
%---------------------------------------------------------------------------------------------------
    function figure_WindowButtonDownFcn(~, ~)
        
        % Resets all the axes titles whenever the window is clicked
        allAxes = findobj(baseTabGroup, 'Type', 'axes');
        for iAx = 1:length(allAxes)
           allAxes(iAx).Title.String = erase(allAxes(iAx).Title.String, ' [SELECTED]'); 
        end
        selected = 0;
    end
%---------------------------------------------------------------------------------------------------
    function image_ButtonDownFcn(src, ~)
        % Append [SELECTED] to the title of any clicked image
        src.Parent.Title.String = [src.Parent.Title.String(end-2:end), ' [SELECTED]'];
        selected = 1;
    end
%---------------------------------------------------------------------------------------------------
    function drawROIButton_Callback(~, ~)
        
        % ---------- Draw ROI and save relevant information about it ----------        
        if selected
            cm = repmat([ rgb('Blue'); rgb('Green'); rgb('Red'); rgb('Cyan'); rgb('Purple'); rgb('Brown'); ...
                rgb('Indigo'); rgb('DarkRed') ; rgb('Magenta') ; rgb('Gold')], 5, 1); % supports up to 50 ROIs
            currcolor = cm(mod(parentIdxROI,size(cm,1))+1,:); % parentIdxROI starts at 1 when gui initializes
            
            % Prompt user to create a polygon ROI
            h = imfreehand(gca);
            pos = getPosition(h);
            xi = pos(:,1);
            yi = pos(:,2);
            myData.ROIs{parentIdxROI}(indexROI).mask = createMask(h);
            myData.ROIs{parentIdxROI}(indexROI).xi = xi;
            myData.ROIs{parentIdxROI}(indexROI).yi = yi;
            currAxes = gca;
            
            % Plot ROI and a numeric identifier
            ROIplots(end+1) = plot(xi, yi, 'linewidth', 2, 'color', currcolor);
            ROItext(end+1) = text(mean(xi),mean(yi),[num2str(parentIdxROI), '.' num2str(indexROI)],'Color',currcolor, 'FontSize',12);
            
            % Save other useful information about the ROI
            myData.ROIs{parentIdxROI}(indexROI).plane = currAxes.UserData;
            myData.ROIs{parentIdxROI}(indexROI).color = currcolor;
            myData.ROIs{parentIdxROI}(indexROI).refImg = refImages(:, :, myData.ROIs{parentIdxROI}(indexROI).plane);
            
            indexROI = indexROI + 1; % Track total # of ROIs that have been drawn for this parent ROI
        else
            disp('Must click to select a plot before drawing an ROI')
        end
    end
%---------------------------------------------------------------------------------------------------
    function addROIButton_Callback(~, ~)
        
        % Update ROI indices
        parentIdxROI = parentIdxROI + 1;
        indexROI = 1;
    end
%---------------------------------------------------------------------------------------------------
    function clearROIButton_Callback(~, ~)
        % Clear all existing ROIs and plots
        for iROI = 1:length(ROIplots)
            try
                delete(ROIplots(iROI))
                delete(ROItext(iROI))
            catch
            end
        end
        myData.ROIs = [];
        indexROI = 1; % Reset global count of total # of ROIs drawn
        parentIdxROI = 1;
        drawnow()
        disp('ROIs cleared')
    end
%---------------------------------------------------------------------------------------------------
    function ROISaveButton_Callback(~,~)
        
        ROImetadata = myData.ROIs;
        
        % Prompt user for save directory
        saveDir = uigetdir(pathName, 'Select a save directory');
        if saveDir == 0
            % Throw error if user canceled without choosing a directory
            disp('ERROR: you must select a save directory or provide one as an argument');
            return
        end
        
        % Prompt user for file name
        fileName = inputdlg('Please choose a file name', 'Save ROI metadata', 1, {'ROI_metadata'});
        fileName = fileName{:};
        
        % Warn user and offer to cancel save if this video will overwrite an existing file
        overwrite = 1;
        if exist(fullfile(saveDir, [fileName, '.mat']), 'file') ~= 0
            dlgAns = questdlg('Creating this data will overwrite an existing file in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
            if strcmp(dlgAns, 'No')
                overwrite = 0;
                disp('Saving cancelled')
            end
        end
        
        % Save ROI data
        if overwrite
            save(fullfile(saveDir, [fileName, '.mat']), 'ROImetadata', '-v7.3');
            disp('ROI data saved!')
        end
        
    end%function

%---------------------------------------------------------------------------------------------------
    function loadPCAButton_Callback(~,~)
        
        % Prompt user to select a file of PCA data (set default dir to match ref image file)
        [pcaFile, pcaFilePath, ~] = uigetfile('*.mat', 'Select a file of PCA data', pathName);
        if pcaFile == 0
            disp('No file selected...PCA data not loaded')
        else
            % Load PCA data
            pcaVars = load(fullfile(pcaFilePath, pcaFile)); % --> 'pcaData' = [x, y, pc, plane], 'explained'
            pcaData = pcaVars.pcaData;
            
            % Enable slider
            slider = findobj('Tag', 'planeNumSlider');
            slider.Enable = 'on';
            slider.Value = 1;
            
            % Call slider callback to plot PCA data for first plane
            planeNumSlider_Callback(findobj('Tag', 'planeNumSlider'), [])
        end
    end

%---------------------------------------------------------------------------------------------------
    function planeNumSlider_Callback(src,~)
        
        % -------- Plot reference image and first 8 PCs for plane #1 -------
        
        % Make sure slider is on an integer value
        src.Value = round(src.Value);
        
        % Calculate necessary axis dimensions
        subplotCols = 2;
        subplotRows = 2;
        padSize = 0.025;
        rightMargin = 0.12;
        axWidth = (1 - rightMargin - ((subplotCols + 1) * padSize)) / subplotCols;
        axHeight = (1 - ((subplotRows + 1) * padSize)) / subplotRows;
        
        % Delete existing pca Axes if necessary
        oldAx = findobj(src.Parent, 'Tag', 'pcaAxes');
        for iAx = 1:numel(oldAx)
            delete(oldAx(iAx));
        end
        
        % Place axes in correct locations
        pcaScreeningAxes = [];
        for yDim = 1:subplotRows
            for xDim = 1:subplotCols
                xPos = (xDim * padSize) + ((xDim-1) * axWidth);
                yPos = (yDim * padSize) + ((yDim-1) * axHeight);
                pcaScreeningAxes{yDim, xDim} = axes(pcaTab, 'Units', 'Normalized', 'Position', ...
                    [xPos, yPos, axWidth, axHeight], 'Tag', 'pcaAxes', 'NextPlot', 'add');
            end
        end
        
        % Order handles by z-planes, displayed from L-R, top-bottom
        pcaScreeningAxes = reshape(flipud(pcaScreeningAxes)', 1, subplotCols * subplotRows);
        
        % Plot reference image in first plot and PCA data in the rest
        pcaImages = [];
        nPlots = subplotCols * subplotRows;
        currPlane = src.Value;
        for iPlot = 1:nPlots
            currAxes = pcaScreeningAxes{iPlot};
            axes(currAxes); axis off; cla(currAxes); currAxes.Title.String = '';
            if iPlot == 1
               pcaRefImg = imshow(refImages(:, :, currPlane), [0 MAX_INTENSITY]);
               currAxes.Title.String = ['Plane #' num2str(src.Value)];
            else
                pcaImages{iPlot - 1} = imagesc((imgaussfilt(pcaData(:,:, iPlot-1, currPlane),0.5)), ...
                                                'Parent', currAxes);
                currAxes.YDir = 'reverse'; % Better than flipud on the input data because that would flip the ROI coordinates as well.
                axis off
                colormap(currAxes, 'bluewhitered')
                pcaImages{iPlot - 1}.ButtonDownFcn = {@image_ButtonDownFcn};
                currAxes.Title.String = ['PC #' num2str(iPlot) - 1];  
                currAxes.UserData = currPlane;
            end
            
        end        
    end


end