
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2019_01_14_exp_2';
sid = 0;

tifFiles = dir(fullfile(parentDir, ['*sid_', num2str(sid), '*.tif']));
controlMeans = []; stimMeans = []; controlDiffs = []; stimDiffs = []; volAvgFrames = [];
for iFile = 1:numel(tifFiles)
    
    % Load data
    tifPath = fullfile(parentDir, tifFiles(iFile).name);
    output = squeeze(read_tif(tifPath));
    disp([num2str(iFile), ' of ', num2str(numel(tifFiles))]);
    frameCounts(iFile) = size(output, 3);
    
    % Divide in half
    nLines = size(output, 1);
    stimROI = output(1:nLines/2, :, :);
    controlROI = output(nLines/2 + 1:end, :, :);
    
    % Get mean values for each frame
    stimMean = squeeze(mean(mean(stimROI, 1), 2));
    controlMean = squeeze(mean(mean(controlROI, 1), 2));
    stimDiff = [0; diff(stimMean)];
    controlDiff = [0; diff(controlMean)];
    
    % Get volume-averaged image
    volAvgFrames(:,:,iFile) = mean(output, 3);
    
    % Save stim frame data
    stimMeans{iFile} = stimMean;
    controlMeans{iFile} = controlMean;
    stimDiffs{iFile} = stimDiff;
    controlDiffs{iFile} = controlDiff;
end


refImg = mean(volAvgFrames, 3);
figure; imshow(refImg, []);

% 
% stimOnFrames = []; stimOffFrames = [];
% % Find stim start and end frames
% [~, stimOnFrame] = max(stimDiff);
% [~, stimOffFrame] = min(stimDiff);
% stimOnFrames(iFile) = stimOnFrame;
% stimOffFrames(iFile) = stimOffFrame;




