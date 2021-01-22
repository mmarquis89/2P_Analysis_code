function create_average_fluorescence_vid(parentDir, sid, varargin)
%===================================================================================================
% MAKE VIDEO OF THE AVERAGE FLUORESCENCE FOR ALL PLANES AND ALL TRIALS
%
% Loads the imaging data for each trial, and computes the average fluorescence during the trial for
% each plane. Then, plots these average images in a grid and save each figure as a frame of a video.
% Can also handle imaging data with multiple channels
%
% INPUTS:
%
%    parentDir = folder with the imaging data files in it
%
%    sid      = a numeric vector specifying the session ID that you would like to process
%
% OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%
%   'OutputDir' = (default: parentDir) the directory to save the output video in
%
%===================================================================================================
try
    addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster
    addpath(genpath('/home/mjm60/DownloadedFunctions')) % if running on O2 cluster
    
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'OutputDir', parentDir);
    parse(p, varargin{:});
    outputDir = p.Results.OutputDir;
    
    % Identify data files for each session
    myFiles = dir(fullfile(parentDir, '*sid*bid*.tif'));
    fileNames = sort({myFiles.name})';
    
    write_to_log(['Found ', num2str(numel(fileNames)), ' data files'], mfilename);
    
    sidLocs = regexp(fileNames, repmat({'(?<=sid_).'}, length(fileNames), 1));
    for iFile = 1:length(fileNames)
        sessionNums(iFile) = str2double(fileNames{iFile}(sidLocs{iFile}));
    end
    currSidFiles = fileNames(sessionNums == sid);
    for iFile = 1:length(currSidFiles)
        
        % Load a .tif file
%         write_to_log(['About to read .tif file: ', fullfile(parentDir, currSidFiles{iFile})], mfilename)
        rawFile = read_tif(fullfile(parentDir, currSidFiles{iFile}));
        rawFile = squeeze(rawFile); % [Lines Pixels Planes Volumes Channels]
%         write_to_log(['Raw tif size = ', num2str(size(rawFile))], mfilename);
        chanData = zeros(size(rawFile)); chanData = chanData(:,:,1:end-4,:,:);
        nChannels = size(chanData, 5);
        nPlanes = size(chanData, 3);
        
        if iFile == 1
            write_to_log('Creating video writer(s)...', mfilename)
            for iChannel = 1:nChannels
                
                % Create save directory if it does not already exist
                if ~isdir(outputDir)
                    mkdir(outputDir)
                end
                
                % Create a video writer for each channel
                fileName = ['Channel_', num2str(iChannel), '_sid_', num2str(sid), '_average_fluorescence_by_trial'];
                myVids(iChannel) = VideoWriter(fullfile(outputDir, fileName));
                myVids(iChannel).FrameRate = 1;
                open(myVids(iChannel));
            end
        end
        
        for iChannel = 1:nChannels
            
            % Offset so minimum value = 1
            workingFile = squeeze(rawFile(:,:,:,:,iChannel));
            workingFile = workingFile-min(workingFile(:));
            workingFile = workingFile + 1;
            
            % Discard flyback frames
            workingFile(:,:,end-(4-1):end,:) = [];                         % --> [y, x, plane, volume]
            chanData(:,:,:,:,iChannel) = workingFile(:,:,:,:,1); % --> [y, x, plane, volume, channel]
%             write_to_log('Flyback frames discarded', mfilename);
        end
        
        
        % Write averaged frames to video
        for iChannel = 1:nChannels
            
            % Average across all volumes
            avgData = squeeze(mean(chanData(:,:,:,:,iChannel), 4));       % --> [y, x, plane]
%             write_to_log('Volume average calculated', mfilename);
            
            % Create fig
            f = figure(1);clf
            f.Position = [50 45, 1620, 950];
            
            % Figure out how many subplots are needed
            nPlots = numSubplots(nPlanes);
            
            for iPlane = nPlanes:-1:1 % Reverse order so planes go from dorsal --> ventral
                
                % Plot averaged image for each plane
                subaxis(nPlots(1), nPlots(2), iPlane, 'Spacing', 0, 'MB', 0.025);
                imshow(squeeze(avgData(:,:,iPlane)), [])
                
                % Label postions
                if iPlane == nPlanes
                    title('Ventral')
                elseif iPlane == 1
                    title('Dorsal')
                end
                
            end%iPlane
            
            % Write frame to video
            drawnow
            writeFrame = getframe(f);
            writeVideo(myVids(iChannel), writeFrame);
            close(f);
            
        end%iChannel
    end%iFile
    
    for iChannel = 1:nChannels
        close(myVids(iChannel))
    end
catch ME
    write_to_log(getReport(ME), mfilename);
end
end