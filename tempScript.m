parentDir = 'D:\Dropbox (HMS)\2P Data\20191107-1_38A11-Chrimson_60D05-7f';
outputDir = fullfile(parentDir, 'ProcessedData');

[test, tifMetadata] = read_tif(fullfile(parentDir, 'Stack_00022.tif'));

tifMetadata = parse_scanimage_metadata(tifMetadata);


% Make anatomy stack
create_anatomy_stack(parentDir, 'FileString', 'Stack_*.tif', 'OutputFilePrefix', 'AnatomyStack', ...
        'OutputDir', outputDir);
    
    
%% CONSOLIDATE METADATA AND DAQ DATA FOR ALL TRIALS

% Process metadata files
mdFiles = dir(fullfile(parentDir, '*metadata*trial*.mat'));
expMetadata = [];
for iFile = 1:numel(mdFiles)
       
    % Load the file 
    % Contains struct "mD" with fields [trialSettings], [expID], [expDirName], [expDir], 
    % [trialNum], and [SAMPLING_RATE])
    load(fullfile(parentDir, mdFiles(iFile).name)); 
        
    % Flatten structure so the contents of trialSettings are at the root level
    mD = rmfield(setstructfields(mD, mD.trialSettings), 'trialSettings');

    % Concatenate into master structure
    if iFile == 1
        expMetadata = mD;
    else
        expMetadata(iFile) = mD;
    end

end%iFile

% Process DAQ data files
daqDataFiles = dir(fullfile(parentDir, '*daqData*trial*.mat'));
expDaqData = [];
for iFile = 1:numel(daqDataFiles)
       
    % Load the file 
    % Contains variabls "outputData", "trialData", and "ColumnLabels"
    currData = load(fullfile(parentDir, daqDataFiles(iFile).name)); 
    
    % Get the trial number of the current file
    currData.trialNum = str2double(regexp(daqDataFiles(iFile).name, '(?<=trial_)...(?=.mat)', ...
            'match', 'once'));
        
    % Add to master structure  
    expDaqData(iFile) = currData;
    
end%iFile

% Save in processed data directory
save(fullfile(outputDir, 'metadata.mat'), 'expMetadata', '-v7.3');
save(fullfile(outputDir, 'daqData.mat'), 'expDaqData', '-v7.3');

%% PROCESS FICTRAC DATA




%% PROCESS IMAGING DATA AND METADATA

% Identify imaging data files
imgDataFiles = dir(fullfile(parentDir, '*trial*.tif'));
for iFile = 1:numel(imgDataFiles)

    disp(['Processing file #', num2str(iFile), ' of ', num2str(numel(imgDataFiles))]);
    
    % Get the trial number of the current file
    trialNum = str2double(regexp(imgDataFiles(iFile).name, '(?<=trial_)...(?=.*\.tif)', ...
            'match', 'once'));
        
    % Load the current file
    [imgData, tifMetadata] = read_tif(fullfile(parentDir, imgDataFiles(iFile).name));
    
    % Parse ScanImage metadata
    siMetadata = parse_scanimage_metadata(tifMetadata);
    siMetadata.trialNum = trialNum;
    if iFile == 1
        imagingMetadata = siMetadata;
    else
        imagingMetadata(iFile) = siMetadata;
    end
    
    % Discard flyback frames from imaging data
    nFlybackFrames = siMetadata.SI.hFastZ.numDiscardFlybackFrames;
    imgData(:, :, (end - (nFlybackFrames - 1)):end, :) = []; % --> [y, x, plane, volume]
    
    % Doing these next two steps because for some reason ScanImage offsets fluorescence data values 
    % by some random amount in different individual .tif files...I've found that this approach fixes 
    % that issue and makes the data directly comparable across trials
    
        % Clip bottom 5% of minimum pixel values per frame
        minFrameVals = sort(as_vector(min(min(imgData))));
        clipThresh = minFrameVals(round(length(minFrameVals) * 0.05));
        imgData(imgData < clipThresh) = clipThresh;

        % Then offset so min value = 1 
        imgData = imgData - min(imgData(:));
        imgData = imgData + 1;    
        
    % Save the raw normalized data in a .mat file
    disp('Saving processed data...')
    save(fullfile(outputDir, ['imagingData_trial_', pad(num2str(trialNum), 3, 'left', '0'), '.mat']), ...
            'imgData', '-v7.3');
            
    % Create and save a file containing reference images for each file
    refImages = mean(imgData, 4); % --> [y, x, plane]
    save(fullfile(outputDir, ['refImages_raw_trial_' pad(num2str(trialNum), 3, 'left', '0'), '.mat']), ...
            'refImages', '-v7.3');
    
end%iFile

% Save the imaging metadata
save(fullfile(outputDir, 'imagingMetadata.mat'), 'imagingMetadata', '-v7.3');



