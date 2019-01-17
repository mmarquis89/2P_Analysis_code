
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2019_01_15_exp_3';
sid = 0;

tifFiles = dir(fullfile(parentDir, ['*sid_', num2str(sid), '*.tif']));
controlMeans = []; stimMeans = []; controlDiffs = []; stimDiffs = []; stimOnFrames = []; stimOffFrames = [];
for iFile = 1:numel(tifFiles)
    
    % Load data
    tifPath = fullfile(parentDir, tifFiles(iFile).name);
    output = squeeze(read_tif(tifPath));
    disp([num2str(iFile), ' of ', num2str(numel(tifFiles))]);
    
    % Divide in half
    nLines = size(output, 1);
    stimROI = output(1:nLines/2, :, :);
    controlROI = output(nLines/2 + 1:end, :, :);
    
    % Get mean values for each frame
    stimMean = squeeze(mean(mean(stimROI, 1), 2));
    controlMean = squeeze(mean(mean(controlROI, 1), 2));
    stimDiff = [0; diff(stimMean)];
    controlDiff = [0; diff(controlMean)];
    
    % Find stim start and end frames
    [~, stimOnFrame] = max(stimDiff);
    [~, stimOffFrame] = min(stimDiff);
    
    % Save a frame from the control and stim period for each trial
    stimFrame{iFile} = output(:,:, stimOnFrame + 1);
    controlFrame{iFile} = output(:,:, stimOffFrame + 1);
    
    % Save stim frame data
    stimOnFrames(iFile) = stimOnFrame;
    stimOffFrames(iFile) = stimOffFrame;
    stimMeans{iFile} = stimMean;
    controlMeans{iFile} = controlMean;
    stimDiffs{iFile} = stimDiff;
    controlDiffs{iFile} = controlDiff;
end

test = [];
for iTrial = 1:numel(stimMeans)
    test(:,iTrial) = stimMeans{iTrial};
end


