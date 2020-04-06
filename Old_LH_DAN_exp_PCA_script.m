
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017\2017_09_27_exp_2';

%%
% PROCESS IMAGING DATA AND METADATA

% Identify imaging data files
imgDataFiles = dir(fullfile(parentDir, 'cdata*.tif'));

for iTrial = 1:numel(imgDataFiles)

    disp(['Processing file #', num2str(iTrial), ' of ', num2str(numel(imgDataFiles))]);
    
    % Get the trial number of the current file
    trialNum = str2double(regexp(imgDataFiles(iTrial).name, '(?<=\__).....(?=.*\.tif)', ...
            'match', 'once'));
        
    % Load the current file
    [imgData, tifMetadata] = read_tif(fullfile(parentDir, imgDataFiles(iTrial).name));
    
    % Parse ScanImage metadata
    siMetadata = parse_scanimage_metadata(tifMetadata);
    siMetadata.trialNum = trialNum;
    if iTrial == 1
        imagingMetadata = siMetadata;
    else
        imagingMetadata(iTrial) = siMetadata;
    end
    
    % Detect number of channels
    if numel(size(imgData)) == 5 
        nChannels = 2;
    else
        nChannels = 1;
    end
    
    % Process only one channel
    fileSuffix = '';
    if nChannels == 2
        imgData = imgData(:, :, :, :, 1);
%         imgData = imgData(:, :, :, :, 2);
%         fileSuffix = '_cyRFP';
    end
    
    % Discard flyback frames from imaging data
    if siMetadata.SI.hFastZ.discardFlybackFrames
        nFlybackFrames = siMetadata.SI.hFastZ.numDiscardFlybackFrames;
        imgData(:, :, (end - (nFlybackFrames - 1)):end, :) = []; % --> [y, x, plane, volume]
    elseif size(imgData, 3) == 16 
        % All these are old experiments with 4 flyback frames, and acquired with an older version
        % of ScanImage that stored the flyback frames at the beginning of the stack instead of the 
        % end.
        imgData(:, :, 1:4, :) = [];
    end
    
    % Doing these next two steps because for some reason ScanImage offsets fluorescence data values 
    % by some random amount in different individual .tif files...I've found empirically that this 
    % approach fixes that issue and makes the data directly comparable across trials
            
%         % Clip bottom 5% of minimum pixel values per frame
%         minFrameVals = sort(as_vector(min(min(imgData))));
%         clipThresh = minFrameVals(round(length(minFrameVals) * 0.05));
%         imgData(imgData < clipThresh) = clipThresh;
% 
%         % Then offset so min value = 1 
%         imgData = imgData - min(imgData(:));
%         imgData = imgData + 1;    
        
    % Save the raw data in a .mat file
    disp('Saving processed data...')
    save(fullfile(parentDir, ['imagingData_raw_trial_', pad(num2str(trialNum), 3, 'left', '0'), ...
            fileSuffix, '.mat']), 'imgData', '-v7.3');
            
%     % Create and save a file containing reference images for each file
%     refImages = mean(imgData, 4); % --> [y, x, plane]
%     save(fullfile(outputDir, ['refImages_raw_trial_' pad(num2str(trialNum), 3, 'left', '0'), ...
%             fileSuffix, '.mat']), 'refImages', '-v7.3');
end

%% Combine all trials into a single array, then save data for individual planes

dataFiles = dir(fullfile(parentDir, 'imagingData_raw_trial*.mat'));
for iFile = 1:numel(dataFiles)
   
    disp(iFile);
    load(fullfile(parentDir, dataFiles(iFile).name), 'imgData');
    
    if iFile == 1
        wholeSession = zeros([size(imgData), numel(dataFiles)]);
    end
    
    wholeSession(:, :, :, :, iFile) = imgData;
    
end
sz = size(wholeSession);
minVal = min(wholeSession(:));

%
for iPlane = 1:sz(3)
    disp(iPlane);
    imgData = squeeze(wholeSession(:, :, iPlane, :, :)) - minVal;
    saveFileName = ['imagingData_raw_plane_', pad(num2str(iPlane), 3, 'left', '0'), '.mat'];
    save(fullfile(parentDir, saveFileName), 'imgData', '-v7.3')
end

clear wholeSession;

%
fileSuffix = '';

% ----- Plane-wise (2D) NoRMCorre motion correction -----

imgDataFiles = dir(fullfile(parentDir, 'imagingData_raw*plane*.mat'));

for iPlane = 1:numel(imgDataFiles)
    
    load(fullfile(parentDir, imgDataFiles(iPlane).name), 'imgData'); % --> [y, x, volume, trial]
    [nLines, nCols, nVolumes, nTrials] = size(imgData);
    
    currData = reshape(imgData, [nLines, nCols, nVolumes*nTrials]); % --> [y, x, trialVolume]
    
    % Clip highest 0.1% of values and smooth with 2D gaussian filter
    srt = sort(currData(:));
    capVal = srt(numel(srt) - round(numel(srt)/1000));
    currData(currData > capVal) = capVal;
    currDataSm = zeros(size(currData));
    for iVol = 1:size(currData, 3)
%         disp(iVol)
        currDataSm(:, :, iVol) = imgaussfilt(currData(:, :, iVol), 0.5);
    end
    currData = currDataSm;
    
    % Set NoRMCorre options
    options_rigid = NoRMCorreSetParms('d1', size(currData, 1), 'd2', size(currData, 2), ...
            'max_shift', [25, 25], ...
            'init_batch', 100, ...
            'us_fac', 50 ...
            );
    
    % Run registration
    tic; [planeData, ~, regTemplate, ~] = normcorre_batch(currData, options_rigid); toc
    
    
    % Save registered data reference images
    refImages = mean(planeData, 3);    % --> [y, x, plane]
    save(fullfile(parentDir, ['refImages_reg_plane_', pad(num2str(iPlane), 3, 'left', '0'), ...
        fileSuffix, '.mat']), 'refImages', 'regTemplate', '-v7.3');
    
    % Save registered imaging data
    imgData = reshape(planeData, [nLines, nCols, nVolumes, nTrials]);
    disp('Saving registered data...');
    save(fullfile(parentDir, ['imagingData_reg_plane_', pad(num2str(iPlane), 3, 'left', '0'), ...
        fileSuffix, '.mat']), 'imgData', '-v7.3')
    

    
end%iPlane

%% Load all registered imaging data and save in blocks of trials

imgDataFiles = dir(fullfile(parentDir, 'imagingData_reg_plane*.mat'));
load(fullfile(parentDir, imgDataFiles(1).name), 'imgData');
sz = size(imgData);

wholeSession = zeros([sz(1:2), numel(imgDataFiles), sz(3:4)]);
[nLines, nCols, nPlanes, nVolumes, nTrials] = size(wholeSession);

for iFile = 1:numel(imgDataFiles)
    disp(iFile)
    load(fullfile(parentDir, imgDataFiles(iFile).name), 'imgData');
    wholeSession(:, :, iFile, :, :) = imgData;
end

filesPerBlock = 20;
nBlocks = sz(4) / filesPerBlock;
for iBlock = 1:nBlocks
    currBlockStart = 1 + (filesPerBlock * (iBlock-1));
    currBlockEnd = currBlockStart + (filesPerBlock - 1);
    imgData = wholeSession(:, :, :, :, currBlockStart:currBlockEnd);
     
    save(fullfile(parentDir, ['imagingData_reg_block_', pad(num2str(iBlock), 3, 'left', '0'), ...
            '.mat']), 'imgData', '-v7.3');
    
    
end


%% Run whole-volume PCA on blocks of trials

blockDataFiles = dir(fullfile(parentDir, 'imagingData_reg_block_*.mat'));

nPCs = 30;
sigma = 0.5;

for iFile = 1:numel(blockDataFiles)
   
    load(fullfile(parentDir, blockDataFiles(iFile).name), 'imgData');
    [nLines, nCols, nPlanes, nVolumes, nTrials] = size(imgData);
    
    imgData = reshape(imgData, [nLines, nCols, nPlanes, (nVolumes * nTrials)]);
    nBlockVolumes = size(imgData, 4);
    
    % Pre-smooth with gaussian filter
    for iVol = 1:nBlockVolumes
        for iPlane = 1:nPlanes
            imgData(:, :, iPlane, iVol) = imgaussfilt(imgData(:, :, iPlane, iVol), sigma);
        end
    end
    sz = size(imgData);
    
    % Run PCA
    data2D = double(reshape(imgData, [(sz(1)*sz(2)*sz(3)), sz(4)]));     % --> [pixel, trialVolume]
    tData2D = data2D';                                                   % --> [trialVolume, pixel]
    tic
    [coeff,score,latent,~,explained] = pca(tData2D, 'numcomponents', nPCs); % [pixel, pc]
    tc = toc;
    pcaData = reshape(coeff, [sz(1:3), nPCs]);                          % [y, x, plane, pc]
    
%     fileName = ['pcaData_trial_', pad(num2str(iFile), 3, 'left', '0'), '.mat'];
    fileName = ['pcaData_block_', pad(num2str(iFile), 3, 'left', '0'), '.mat'];
    save(fullfile(parentDir, fileName), 'pcaData')
    disp(['Block ', num2str(iFile), ' PCA completed in ', num2str(tc), ' sec'])
    
    
    
end

%%


blockNum = 4


nPCs = 9;
offset = 0;
targetPCs = [];
targetPlanes = [2:9];
sigma = 0.5;
maxIntensity = 1000;

% Load PCA data for target block
load(fullfile(parentDir, ['pcaData_block_', pad(num2str(blockNum), 3, 'left', '0'), '.mat']), ...
        'pcaData')


% Load full-experiment reference images for each plane
allRefImages = [];
for iPlane = 1:size(pcaData, 3)
    load(fullfile(parentDir, ['refImages_reg_plane_', pad(num2str(iPlane), 3, 'left', '0'), ...
            '.mat']), 'refImages');
    allRefImages(:, :, iPlane) = refImages;
end
refImages = allRefImages;

if ~isempty(targetPCs)
    nPCs = numel(targetPCs);
end
if ~isempty(targetPlanes)
    nPlanes = numel(targetPlanes);
else
    nPlanes = size(pcaData, 3);
    targetPlanes = 1:nPlanes;
end
f = figure(2); clf;
f.Color = [1 1 1];
subaxis(nPCs + 1, nPlanes, 1, 1, ...
        'm', 0.01, ...
        'p', 0, ...
        'sh', 0.001, 'sv', 0)
    
for iPlane = 1:nPlanes
    
    % Plot reference images
    subaxis(nPCs + 1, nPlanes, iPlane, 1)
    imshow(refImages(:, :, targetPlanes(iPlane)), [0 2000]);
    colormap(gca, 'gray');
    
    % Plot top PCs
    for iPC = 1:nPCs
       
       if isempty(targetPCs)
           currData = squeeze(pcaData(:, :, targetPlanes(iPlane), iPC + offset));
       else
           currData = squeeze(pcaData(:, :, targetPlanes(iPlane), targetPCs(iPC)));
       end
       currData = flipud(currData);
       
       % Plot reference image underneath
       ax1 = subaxis(nPCs + 1, nPlanes, iPlane, iPC + 1); hold on
       h1 = imagesc(ax1, flipud(refImages(:, :, targetPlanes(iPlane)))); axis image
       ax1.CLim = [0 maxIntensity];
       colormap(ax1, 'gray')
       
       ax2 = axes; hold on;
       ax2.Position = ax1.Position;
       
       currData = imgaussfilt(currData, sigma);
       
       h2 = imagesc(ax2, currData, 'alphadata', ~isnan(currData));
       colormap(ax2, 'bluewhitered'); axis image
       ax2.CLim = ax2.CLim * 0.7;
       
       ax1.Visible = 'off';
       ax2.Visible = 'off';
       
    end
    

end


    
    
    
    
    
    