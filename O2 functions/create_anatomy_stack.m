function create_anatomy_stack(dirPath, varargin)
%===================================================================================================
% MAKE AVERAGED ANATOMY STACK FROM INDIVIDUAL 2P VOLUMES
% Will average together all stacks in the directory that meet the requirement specified by 'fileStr'
% and save as .tif files the mean stack as well as a max Z-projection.
%
% INPUTS:
%    dirPath = folder with stack files, e.g. 'B:\Dropbox (HMS)\2P Data\Imaging Data\2017_09_05'
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%    FileString       = (default: '*Stack_*.tif') filter string for dir() selection to get anatomy
%                       stack files.
%
%    OutputFilePrefix = (default: '') name to prepend to the names of output files.
%
%    OutputDir        = (default: dirPath) directory to save the output files in
%
%===================================================================================================
    try
        addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster
    catch
    end
    
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'FileString', '*Stack_*.tif');
    addParameter(p, 'OutputFilePrefix', '');
    addParameter(p, 'OutputDir', dirPath);
    parse(p, varargin{:});
    fileStr = p.Results.FileString;
    outputFilePrefix = p.Results.OutputFilePrefix;
    outputDir = p.Results.OutputDir;
    
    % Get stack files
    stacks = dir(fullfile(dirPath, fileStr));
    nStacks = length(stacks);
    
    % Add underscore to prefix if one is provided
    if ~isempty(outputFilePrefix)
        outputFilePrefix = [outputFilePrefix, '_'];
    end
    
    % Get sum of all stacks
    summedStacks = [];
    for iVol = 1:nStacks
        disp(['Processing stack #', num2str(iVol), ' of ', num2str(nStacks)]);
        if iVol == 1
            firstStack = int16(read_tif(fullfile(dirPath,stacks(iVol).name)));                  % --> [y, x, plane, volume, channel]
            summedStacks = int16(firstStack);                                                   % --> [y, x, plane, volume, channel]
        else
            summedStacks = summedStacks + int16(read_tif(fullfile(dirPath,stacks(iVol).name))); % --> [y, x, plane, volume, channel]
        end
    end

    % Convert to double
    summedStacks = squeeze(double(summedStacks));                                                % --> [y, x, plane, channel]
    
    % Check whether data has multiple channels
    nChannels = size(summedStacks, 4);
    
    % Create save directory if it does not already exist
    if ~isdir(outputDir)
        mkdir(outputDir)
    end


    if nChannels == 1
        
        % Calculate mean by dividing by total number of stacks and scale from 0-1
        avgStack = summedStacks ./ nStacks;             % --> [y, x, plane]
        avgStack = avgStack - min(avgStack(:));         % --> [y, x, plane]
        avgStackScaled = avgStack ./ max(avgStack(:));  % --> [y, x, plane]
        
        % Get max Z-projection of averaged stack and scale intensity range to 0-1
        maxZ = squeeze(max(avgStack, [], 3));           % --> [y, x]
        maxZScaled = maxZ ./ max(maxZ(:));              % --> [y, x]
        
        % Write averaged stack and Z-projection to .tif files
        imwrite(maxZScaled, fullfile(outputDir, [outputFilePrefix, 'MeanMaxZ.tif']));
        imwrite(avgStackScaled(:,:,1), fullfile(outputDir, [outputFilePrefix, 'MeanStack.tif']));
        for iPlane = 2:size(avgStackScaled, 3)
            imwrite(avgStackScaled(:,:,iPlane), fullfile(outputDir, [outputFilePrefix, 'MeanStack.tif']), 'writemode','append');
        end
        
    else % nChannels == 2
        
        % Calculate mean by dividing by total number of stacks
        avgStack = summedStacks ./ nStacks;                  % --> [y, x, plane, channel]
        
        % Separate data by channel and scale intensity range to 0-1
        avgStack_1 = avgStack(:,:,:,1);                      % --> [y, x, plane]
        avgStack_2 = avgStack(:,:,:,2);                      % --> [y, x, plane]
        avgStack_1 = avgStack_1 - min(avgStack_1(:));        % --> [y, x, plane]
        avgStack_2 = avgStack_2 - min(avgStack_2(:));        % --> [y, x, plane]
        avgStackScaled_1 = avgStack_1 ./ max(avgStack_1(:)); % --> [y, x, plane]
        avgStackScaled_2 = avgStack_2 ./ max(avgStack_2(:)); % --> [y, x, plane]
        
        % Get max Z-projection of averaged stack and scale intensity range to 0-1
        maxZ_1 = squeeze(max(avgStack_1, [], 3));            % --> [y, x]
        maxZ_2 = squeeze(max(avgStack_2, [], 3));            % --> [y, x]
        maxZScaled_1 = maxZ_1 ./ max(maxZ_1(:));             % --> [y, x]
        maxZScaled_2 = maxZ_2 ./ max(maxZ_2(:));             % --> [y, x]
        
        % Write Z-projection to .tif files
        imwrite(maxZScaled_1, fullfile(outputDir, [outputFilePrefix, 'MeanMaxZ_1.tif']));
        imwrite(maxZScaled_2, fullfile(outputDir, [outputFilePrefix, 'MeanMaxZ_2.tif']));
        
        % Write averaged stack to .tif files
        imwrite(avgStackScaled_1(:,:,1), fullfile(outputDir, [outputFilePrefix, 'MeanStack_1.tif']));
        for iPlane = 2:size(avgStackScaled_1, 3)
            imwrite(avgStackScaled_1(:,:,iPlane), fullfile(outputDir, [outputFilePrefix, 'MeanStack_1.tif']), 'writemode','append');
        end
        imwrite(avgStackScaled_2(:,:,1), fullfile(outputDir, [outputFilePrefix, 'MeanStack_2.tif']));
        for iPlane = 2:size(avgStackScaled_2, 3)
            imwrite(avgStackScaled_2(:,:,iPlane), fullfile(outputDir, [outputFilePrefix, 'MeanStack_2.tif']), 'writemode','append');
        end
        
    end%if