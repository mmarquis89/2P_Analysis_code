function preReg_routine_MM(parentDir, sid, expDate, varargin)
%===================================================================================================
% Load and process raw ScanImage imaging data
%
% INPUTS:
%       parentDir   = directory containing the raw data that you want to process.
%
%       sid         = the session ID you want to process.
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%
%       'OutputDir' = (default: parentDir) the directory to save the processed data in
%
% NOTE: the files should be sorted in chronological order due to the timestamp at the beginning
% of the filename. If filenames do not sort in this pattern, they must be renamed before processing.
%===================================================================================================
try
    addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster
    
    write_to_log('Starting pre-registration processing...', mfilename);
    
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'OutputDir', parentDir);
    parse(p, varargin{:});
    outputDir = p.Results.OutputDir;
    
    write_to_log('Identifying data files', mfilename);
    
    % Identify data files for each session
    myFiles = dir(fullfile(parentDir, '*sid*bid*.tif'));
    fileNames = sort({myFiles.name})';
    sidLocs = strfind(fileNames, 'sid_');
    for iFile = 1:length(fileNames)
        sessionNums(iFile) = str2double(fileNames{iFile}(sidLocs{iFile}+4));
    end
    
    disp(['Processing ', num2str(numel(fileNames)), ' files...'])
    write_to_log(['Processing ', num2str(numel(fileNames)), ' files...'], mfilename)
    
    % Make session folder for new files if necessary
    if ~isdir(outputDir)
        mkdir(outputDir)
    end
    
    write_to_log('Getting names of files from current session...', mfilename);
    
    % Separate names of files from the current session
    currFiles = fileNames(sessionNums == sid);
    nTrials = numel(currFiles);
    disp(nTrials)
    for iFile = 1:nTrials
        
        write_to_log(['Processing trial #', num2str(iFile)], mfilename);
        disp(['Processing trial #', num2str(iFile)])
        write_to_log(['Attempting to load ', fullfile(parentDir, currFiles{iFile})], mfilename);
        
        % Load the next .tif (or .mat) data file
        rawFile = read_tif(fullfile(parentDir, currFiles{iFile}));
        rawFile = squeeze(rawFile); % [Lines Pixels Planes Volumes Channels]
        nVolumes = size(rawFile, 4);
        nChannels = size(rawFile, 5);

        % Extract scanimage data from the first trial
        if iFile == 1
            write_to_log(exist(fullfile(parentDir, 'cdata_20190909_125626_sid_0_bid_0_dur_1200_nTrials_60_00001.tif'), 'file'), mfilename);
            [~, tifMetadata] = read_patterned_tifdata(fullfile(parentDir, currFiles{iFile}));
            scanimageInfo = tifMetadata.tifinfo;
            nFlybackFrames = frameStringKeyLookup(scanimageInfo.Software, ...
                    'hFastZ.numDiscardFlybackFrames');
            nPlanes = frameStringKeyLookup(scanimageInfo.Software, ...
                    'hFastZ.numFramesPerVolume') - nFlybackFrames;
            nChannels = length(frameStringKeyLookup(scanimageInfo.Software, ...
                    'scanimage.SI.hChannels.channelSave'));

            fclose(tifMetadata.fid);
            
            % Create a .mat file for each plane/channel
            matfiles = [];
            for iChan = 1:nChannels
                for iPlane = 1:nPlanes
                    matfiles{iChan, iPlane} = matfile(fullfile(outputDir, ['sid_', num2str(sid), '_chan_', ...
                            num2str(iChan), '_plane_', num2str(iPlane), '_sessionFile.mat']));
                    matfiles{iChan, iPlane}.sessionData = [];
                end
            end
        end
                
        % Process data and append to each file
        for iChan = 1:nChannels
            
            workingFile = squeeze(rawFile(:,:,:,:,iChan));
            
            % Discard flyback frames
            workingFile(:,:,(end-(nFlybackFrames-1)):end,:) = [];   % --> [y, x, plane, volume]
            
            % Clip bottom 5% of minimum pixel values per frame
            minFrameVals = sort(as_vector(min(min(workingFile))));
            clipThresh = minFrameVals(round(length(minFrameVals) * 0.05));
            workingFile(workingFile < clipThresh) = clipThresh;
            
            % Then offset so min value = 1
            workingFile = workingFile - min(workingFile(:));
            workingFile = workingFile + 1;
                        
            % Save data for current plane
            for iPlane = 1:nPlanes
                currPlaneData = uint16(squeeze(workingFile(:,:,iPlane,:, iChan))); % --> [y, x, volume]
                sz = size(matfiles{iChan, iPlane}, 'sessionData');
                if numel(sz) < 3
                    matfiles{iChan, iPlane}.sessionData = currPlaneData;
                else
                    matfiles{iChan, iPlane}.sessionData(:,:, (sz(3) + 1):(sz(3) + nVolumes)) = ...
                        currPlaneData;
                end
            end
        end
               
    end%iFile
    
    % Save other imaging metadata variables in separate file
    save(fullfile(outputDir, 'imgMetadata'), 'expDate', 'scanimageInfo', 'nPlanes', 'nVolumes', ...
            'nTrials', '-v7.3');
       
%     % Create reference images
%     refImages = [];
%         
%         % Use red channel for reference images if available, otherwise GCaMP channel
%         for iPlane = 1:nPlanes
%            currPlaneData = matfiles{nChannels, iPlane}.sessionData; % --> [y, x, volume]
%            refImages(:,:, iPlane) = mean(currPlaneData, 3);         % --> [y, x, volume]
%         end
%         channelNum = nChannels;
%         save(fullfile(outputDir, sprintf('sid_%.0f_refImages.mat', sid)), 'refImages', ...
%                 'channelNum', '-v7.3');
        
%     write_to_log('Reference images created and saved', mfilename)
catch ME
    write_to_log(getReport(ME), mfilename);
end%try
end%function