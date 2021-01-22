
%% INTIAL PROCESSING OF RAW SCANIMAGE FILES

% % Load exp list
% expList = load_expList('groupName', 'gaplessAcq');
% % expList = expList(~cellfun(@isempty, regexp(expList.expID, '2018.*', 'match', 'once')), 1);

expList = table({'20200901-1'; '20200901-1'}, 'variablenames', {'expID'}); 

for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID);
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', currExpID(end)]; 
    parentDir = find_parent_dir(currExpID);
    
    if exist(fullfile(parentDir, expDirName), 'dir')
        cdataFiles = dir(fullfile(parentDir, expDirName,'cdata*.tif'));
        if ~isempty(cdataFiles)
            medVals = {};
            minVals = {};
            minFrameVals = {};
            volCounts = [];
            for iFile = 1:numel(cdataFiles)
                
                if ~mod(iFile, 20)
                    disp(iFile);
                end
                
                % Load current file
                currFile = fullfile(parentDir, expDirName, cdataFiles(iFile).name);
                [imgData, siMetadata] = read_tif(currFile, 0); % --> [y, x, planes, volumes]
                sz = size(imgData);
                volCounts(iFile) = size(imgData, 4);
                
                
                 % Detect number of channels, if there are >1 process only the first channel 
                if numel(sz) == 5
                    imgData = imgData(:, :, :, :, 1);
                end        
                
                % Compare exp date to first Berg-2 (and thus SI 2019) experiment, and then
                % discard the appropriate flyback frames (storage location changed in new
                % ScanImage version)
                expDateNum = str2double(currExpID(1:8));
                if expDateNum < 20181006 
                    oldScope = 1;
                    nFlybackFrames = 4; % Should be true for all Berg-1 experiments, but is incorrect in the pre-ScanImage 2019 metadata for some reason
                    imgData(:, :, 1:nFlybackFrames, :) = []; % --> [y, x, plane, volume]
                else
                    oldScope = 0;
                    siMetadata = parse_scanimage_metadata(siMetadata);
                    nFlybackFrames = siMetadata.SI.hFastZ.numDiscardFlybackFrames;
                    imgData(:, :, (end - (nFlybackFrames - 1)):end, :) = []; % --> [y, x, plane, volume]
                end
                
                % Calculate some volume-wide summary statistics for comparison across trials
                minFrameVals{iFile} = sort(as_vector(min(min(imgData))));
                rsData = reshape(permute(imgData, [4 1 2 3]), size(imgData, 4), []); % --> [volume, px]
                medVals{iFile} = median(rsData, 2); 
                minVals{iFile} = min(rsData, [], 2); 
                
                % Save as a .mat file
                save(regexprep(currFile, '\.tif', '\.mat'), 'imgData');
            end
            
            % Save summary statistics for all files
            save(fullfile(parentDir, expDirName, 'summaryStats.mat'), 'minFrameVals', 'medVals', ...
                    'minVals')
                
            % Save list of volume counts
            save(fullfile(parentDir, expDirName, 'volCounts.mat'), 'volCounts');
            
        else
            disp(['No imaging data files found in ', fullfile(parentDir, expDirName)])
        end
    else
        disp(['Could not find experiment directory ', fullfile(parentDir, expDirName)])
    end  
end%iExp


%% PLANE-WISE MOTION CORRECTION OF 20-TRIAL CHUNKS OF DATA

% % Load exp list
% expList = load_expList('groupName', 'gaplessAcq');
% % expList = expList(~cellfun(@isempty, regexp(expList.expID, '2018.*', 'match', 'once')), 1);


for iExp = 1:size(expList, 1)
    currExpID = expList.expID{iExp};
    disp(currExpID);
    expDirName = [currExpID(1:4), '_', currExpID(5:6), '_', currExpID(7:8), '_exp_', ...
            currExpID(end)]; 
    parentDir = find_parent_dir(currExpID);
    
    if exist(fullfile(parentDir, expDirName), 'dir')
        imgDataFiles = dir(fullfile(parentDir, expDirName,'cdata*.mat'));
        if ~isempty(imgDataFiles)
            fileNames = {imgDataFiles.name}';
            
            % Extract their sid, block numbers, and trial numbers
            sidNums = str2double(regexp(fileNames, '(?<=sid_).', 'match', 'once'));
            bidNums = str2double(regexp(fileNames, '(?<=bid_).', 'match', 'once'));
            trialNums = str2double(regexp(fileNames, '(?<=_).....(?=\.mat)', 'match', 'once'));
            
            trialList = table(sidNums, bidNums, trialNums, fileNames, 'VariableNames', ...
                    {'sid', 'blockNum', 'trialNum', 'fileName'});
            
            % Get list of blocks
            acqBlockList = unique(table(sidNums, bidNums, 'VariableNames', {'sid', 'blockNum'}));
            
            for iAcq = 1:size(acqBlockList, 1)
                currTrialList = innerjoin(acqBlockList(iAcq, :), trialList);
                currFileNames = currTrialList.fileName;
                blockStartTrials = 1:20:numel(currFileNames);
                nBlocks = numel(blockStartTrials);
                for iBlock = 1:nBlocks
                    % Identify files in current block
                    disp([expDirName, ', block ', num2str(iBlock), ' of ', num2str(nBlocks)])
                    if iBlock < nBlocks
                        currBlockTrials = blockStartTrials(iBlock):(blockStartTrials(iBlock + 1) ...
                                - 1);
                    else
                        currBlockTrials = blockStartTrials(iBlock):numel(currFileNames);
                    end
                    
                    % Load and concatenate imgData for current block
                    for iFile = 1:numel(currBlockTrials)
                        fileName = fullfile(parentDir, expDirName, ...
                                currFileNames{currBlockTrials(iFile)});
                        disp(['Loading ', currFileNames{currBlockTrials(iFile)}]);
                        load(fileName, 'imgData');
                        if iFile == 1
                            blockData = imgData;
                        else
                            if size(imgData, 4) > size(blockData, 4)
                                blockData = cat(5, blockData, imgData(:, :, :, 1:end-1, :));
                            elseif size(imgData, 4) < size(blockData, 4)
                                blockData = blockData(:, :, :, 1:end-1, :);
                                blockData = cat(5, blockData, imgData);
                            else
                                blockData = cat(5, blockData, imgData); % --> [y, x, plane, volume, trial]
                            end
                        end
                    end
                    
                    % Register data separately for each plane
                    imgDataReg = zeros(size(blockData));
                    regTemplates = [];
                    for iPlane = 1:size(imgDataReg, 3)
                        
                        currData = squeeze(blockData(:,:, iPlane, :)); % --> [y, x, volume]
                        
                        % Clip highest 0.1% of values and smooth with 2D gaussian filter
                        srt = sort(currData(:));
                        capVal = srt(numel(srt) - round(numel(srt)/1000));
                        currData(currData > capVal) = capVal;
                        currDataSm = zeros(size(currData));
                        for iVol = 1:size(currData, 3)
                            currDataSm(:, :, iVol) = imgaussfilt(currData(:, :, iVol), 0.5);
                        end
                        currData = currDataSm;
                        
                        % Set NoRMCorre options
                        options_rigid = NoRMCorreSetParms('d1', size(currData, 1), 'd2', ...
                            size(currData, 2), ...
                            'max_shift', [25, 25], ...
                            'init_batch', 100, ...
                            'us_fac', 50 ...
                            );
                        
                        % Run registration
                        [planeData, ~, regTemplate, ~] = normcorre_batch(currData, options_rigid);
                        
                        % Add data to full volume arrays
                        imgDataReg(:, :, iPlane, :) = planeData;    % --> [y, x, plane, volume, trial]
                        regTemplates(:, :, iPlane) = regTemplate;   % --> [y, x, plane]
                        
                    end%iPlane
                    
                    imgData = imgDataReg;
                    
                    % Save reference images for registered image data
                    refImages = mean(imgData, 4); % --> [y, x, plane]
                    save(fullfile(parentDir, expDirName, ['refImages_reg_block_', ...
                        pad(num2str(iAcq), 3, 'left', '0'), ...
                        '-', num2str(iBlock), '.mat']), 'refImages', ...
                        'regTemplates');
                    
                    % Convert back to int16 data type to conserve storage space
                    imgData = int16(imgData);
                    
                    % Save registered data for the current block
                    disp('Saving registered data...')
                    saveFileName = ['imagingData_reg_block_', pad(num2str(iAcq), 3, ...
                        'left', '0'), '-', num2str(iBlock), '.mat'];
                    save(fullfile(parentDir, expDirName, saveFileName), 'imgData');
                    
                    
                    
                end%iBlock
            end%iAcq
        else
            disp(['No imaging data files found in ', fullfile(parentDir, expDirName)])
        end
    else
        disp(['Could not find experiment directory ', fullfile(parentDir, expDirName)])
    end  
end%iExp






