
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2019_01_11_exp_1';
sid = 0;

tifFiles = dir(fullfile(parentDir, ['*sid_', num2str(sid), '*.tif']));
controlMeans = []; stimMeans = []; controlDiffs = []; stimDiffs = []; stimOnFrames = []; stimOffFrames = [];
for iFile = 1:numel(tifFiles)
%     tifName = 'cdata_20190111_161850_sid_0_bid_4_dur_800_nTrials_40_00001.tif';
    tifPath = fullfile(parentDir, tifName);
    disp(num2str(iFile));
    % Load data
    output = squeeze(read_tif(tifPath));
    
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
    
    % % Plot stuff
    % figure(1); clf; hold on; plot(stimMean), plot(controlMean)
    % figure(2); clf; hold on; plot(stimDiff), plot(controlDiff);
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


