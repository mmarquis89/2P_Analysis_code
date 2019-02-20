function preReg_routine_MM(parentDir, sid, expDate, varargin)
%===================================================================================================
% Load and process raw ScanImage imaging data
%
% INPUTS:
%       parentDir   = directory containing the raw data that you want to process.
%
%       sid        = the session ID you want to process.
%
%       expDate    = the identifier of the experiment you want to process.
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
    
    disp(['Processing', num2str(numel(fileNames)), ' files...'])
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
        
        % Extract scanimage data from the first trial
        if iFile == 1
            [~, tifMetadata] = read_patterned_tifdata(fullfile(parentDir, currFiles{iFile}));
            scanimageInfo = tifMetadata.tifinfo;
            nFlybackFrames = frameStringKeyLookup(scanimageInfo.Software, 'hFastZ.numDiscardFlybackFrames');
            fclose(tifMetadata.fid);
        end
        
%         write_to_log('Tiff data extracted', mfilename);
        
        % Load the .tif file
        rawFile = read_tif(fullfile(parentDir, currFiles{iFile}));
        rawFile = squeeze(rawFile); % [Lines Pixels Planes Volumes Channels]
        chanData = zeros(size(rawFile)); chanData = chanData(:,:,1:end-nFlybackFrames,:,:);
        nChannels = size(chanData, 5);
        
%         write_to_log(['Raw tiff data loaded...dims = ', num2str(size(chanData))], mfilename);
        
        for iChannel = 1:nChannels
            
            % Offset so minimum value = 1
            workingFile = squeeze(rawFile(:,:,:,:,iChannel));
            workingFile = workingFile - min(workingFile(:));
            workingFile = workingFile + 1;
            
            %             % Crop edges (why do I even do this here?? Need to ask AKM what this part was for)
            %             yrange = 3:size(workingFile,1)-2;
            %             xrange = 3:size(workingFile,2)-2;
            %             croppedFile = workingFile(yrange,xrange,:,:);
            
            % Discard flyback frames
            workingFile(:,:,(end-(nFlybackFrames-1)):end,:) = [];       % --> [y, x, plane, volume]
%             write_to_log('Flyback frames trimmed', mfilename);
%             write_to_log(num2str(size(workingFile)), mfilename);
%             write_to_log(num2str(size(chanData)), mfilename);
            chanData(:,:,:,:,iChannel) = workingFile(:,:,:,:,1);    % --> [y, x, plane, volume, channel]
        end
        
        % Separate channels if file includes data from both PMTs
        if nChannels == 1
            chanData_1 = squeeze(chanData);
        else
            chanData_1 = squeeze(chanData(:,:,:,:,1));
            chanData_2 = squeeze(chanData(:,:,:,:,2));
        end
        
        %----- Save to session structure -----
        
        % Create session array(s) on first loop
        if iFile == 1
            chDataSize = size(chanData_1);
            sz = [chDataSize(1:4), nTrials];
            wholeSession = zeros(sz);
            if nChannels > 1
                wholeSession2 = zeros(sz);
            end
        end
        
        fName = currFiles{iFile};

        % Save data to session array(s)
%         write_to_log(['iFile = ', num2str(iFile), ', size(wholeSession) = ', ...
%                     num2str(size(wholeSession)), ', size(chanData_1) = ', ...
%                     num2str(size(chanData_1))], mfilename);
        wholeSession(:,:,:,:,iFile) = chanData_1(:,:,:,:,1);
        if nChannels > 1
            wholeSession2(:,:,:,:,iFile) = chanData_2(:,:,:,:,2);
        end
        
    end%iFile
    
    disp(size(wholeSession))
    
    write_to_log('Session data processing complete', mfilename)
    
    % Save session data file
    save(fullfile(outputDir, ['sid_', num2str(sid), '_sessionFile.mat']), 'wholeSession', '-v7.3')
    if nChannels == 2
        clear wholeSession
        wholeSession = wholeSession2;
        save(fullfile(outputDir, [sid_', num2str(sid), '_Chan_2_sessionFile.mat']), 'wholeSession', '-v7.3');
    end
    
    disp(outputDir)
    disp('Session data file saved')
    
    % Save other imaging metadata variables in separate file
    save(fullfile(outputDir, 'imgMetadata'), 'expDate', 'scanimageInfo', '-v7.3');
    
    disp('Metadata saved')
    
    write_to_log('Session data saved', mfilename);
    
    % Calculate and save volume-averaged session data
    volAvgSessionData = squeeze(mean(wholeSession, 4)); % --> [y, x, plane, trial]
    save(fullfile(outputDir, sprintf('sid_%.0f_volAvgSessionData.mat', sid)), 'volAvgSessionData', '-v7.3');
    
    % Create reference images
    if nChannels == 2
        
        % Use red channel for reference images
        channelNum = 2;
        refImages = [];
        for iPlane = 1:sz(3)
            refImages{iPlane} = squeeze(mean(mean(wholeSession2(:,:,iPlane,:,:),4),5)); % --> [y, x]
        end
        save(fullfile(outputDir, sprintf('sid_%.0f_refImages.mat', iSession)), 'refImages', 'channelNum', '-v7.3');
        
    else
        
        % Use the GCaMP channel to create and save reference images for each plane and for each trial
        channelNum = 1;
        refImages = [];
        for iPlane = 1:sz(3)
            refImages{iPlane} = squeeze(mean(volAvgSessionData(:,:,iPlane,:),4)); % --> [y, x]
        end
        save(fullfile(outputDir, sprintf('sid_%.0f_refImages.mat', sid)), 'refImages', 'channelNum', '-v7.3');
    end
    
    write_to_log('Reference images created and saved', mfilename)
catch ME
    write_to_log(getReport(ME), mfilename);
end%try
end%function