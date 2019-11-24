
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\20191119-1_38A11_ChR_60D05_7f\ProcessedData';

% Load ROI data
disp('Loading data...')
load(fullfile(parentDir, 'roiData_reg.mat'), 'allROIData');

% Load FicTrac data
load(fullfile(parentDir, 'FicTracData.mat'), 'ftData');

% Identify reference image files
refImgFileStr = 'refImages_reg';
refImgFiles = dir(fullfile(parentDir, [refImgFileStr, '*.mat']));
refImgFileNums = get_trialNum({refImgFiles.name});

% Load metadata files
load(fullfile(parentDir, 'daqData.mat'), 'expDaqData');
load(fullfile(parentDir, 'metadata.mat'), 'expMetadata');
load(fullfile(parentDir, 'imagingMetadata.mat'), 'imagingMetadata');

% COMPILE METADATA

allTrialNums = unique([[allROIData.trialNum], [expDaqData.trialNum], [expMetadata.trialNum], ...
        [imagingMetadata.trialNum], [ftData.trialNum]]); 

mD = []; 
for iTrial = 1:(max(allTrialNums) - min(allTrialNums) + 1)
    
    disp(iTrial);
    
    tid = allTrialNums(iTrial);
    mD(iTrial).trialNum = tid;
    
    % Imaging metadata
    if sum([imagingMetadata.trialNum] == tid)
        mD(iTrial).SI = imagingMetadata([imagingMetadata.trialNum] == tid).SI;
    else
        mD(iTrial).SI = [];
    end
    
    % Recorded panels position data
    if sum([expDaqData.trialNum] == tid)
        mD(iTrial).panelsPosData = expDaqData([expDaqData.trialNum] == tid).trialData(:, 1:2); % --> [sample, dim]
    else
        mD(iTrial).panelsPosData = [];
    end
    
    % ROI data
    if sum([allROIData.trialNum] == tid)
        currROIData = allROIData([allROIData.trialNum] == tid).roiDefs;
        for iROI = 1:numel(currROIData)
            
            % Extract just the most useful components of each ROI
            mD(iTrial).roiData(iROI).name = currROIData(iROI).name;
            mD(iTrial).roiData(iROI).color = currROIData(iROI).color;
            for iSubROI = 1:numel(currROIData(iROI).subROIs)
                currSubROI = currROIData(iROI).subROIs;
                mD(iTrial).roiData(iROI).subROIs(iSubROI).plane = currSubROI.plane;
                mD(iTrial).roiData(iROI).subROIs(iSubROI).position = currSubROI.position;
            end
            mD(iTrial).roiData(iROI).rawData = currROIData(iROI).rawData;
            mD(iTrial).roiData(iROI).dffData = currROIData(iROI).dffData;
            mD(iTrial).roiData(iROI).zscoreData = currROIData(iROI).zscoreData;
        end
        
        % Sort ROIs alphabetically by name
        [~, sortIdx] = sort({mD(iTrial).roiData.name});
        mD(iTrial).roiData = mD(iTrial).roiData(sortIdx);
    else
        mD(iTrial).roiData = [];
    end
    
    % Reference images
    currFile = refImgFiles(refImgFileNums == tid);
    if ~isempty(currFile)
        load(fullfile(parentDir, currFile.name), 'refImages');
        mD(iTrial).refImages = refImages; % --> [y, x, plane]
    else
        mD(iTrial).refImages = [];
    end
    
    % FicTrac data
    if sum([ftData.trialNum] == tid)
        currFtData = ftData([ftData.trialNum] == tid).trialData;
        
        % Get frame times relative to start of trial (column 24 is nanosec since previous frame)
        mD(iTrial).ftFrameTimes = cumsum(currFtData(:, 24)) ./ 1e9;
        
        % Copy important variables, converting units as needed
        IFI = [mD(iTrial).ftFrameTimes(1); diff(mD(iTrial).ftFrameTimes)];
        mD(iTrial).ftData.intX = currFtData(:, 15) * 4.5;               % mm
        mD(iTrial).ftData.intY = currFtData(:, 16) * 4.5;               % mm
        mD(iTrial).ftData.intHD = currFtData(:, 17);                    % radians
        mD(iTrial).ftData.moveSpeed = (currFtData(:, 19) * 4.5) ./ IFI; % mm/sec
        mD(iTrial).ftData.intFwMove = currFtData(:, 20) * 4.5;          % mm
        mD(iTrial).ftData.intSideMove = currFtData(:, 21) * 4.5;        % mm
        
        % Calculate derived FicTrac variables      
        mD(iTrial).ftData.yawSpeed = [0; diff(smoothdata(unwrap(currFtData(:, 17), [], 1), ...
                1, 'gaussian', 7), 1)] ./ IFI;  % radians/sec
        mD(iTrial).ftData.fwSpeed = [0; diff(smoothdata(mD(iTrial).ftData.intFwMove, 1, ...
                'gaussian', 7), 1)] ./ IFI;     % mm/sec
        mD(iTrial).ftData.sideSpeed = [0; diff(smoothdata(mD(iTrial).ftData.intSideMove, 1, ...
                'gaussian', 7), 1)] ./ IFI;     % mm/sec      
            
    else
        mD(iTrial).ftData = [];
    end    
    
    % Experiment parameters
    if sum([expMetadata.trialNum] == tid)
        allFieldNames = fieldnames(expMetadata);
        topFields = {'expID', 'SAMPLING_RATE', 'trialDuration', 'using2P', 'usingPanels', ...
                'displayRate', 'usingOptoStim', 'optoStimTiming'};
        for iField = 1:numel(allFieldNames)
            currField = allFieldNames{iField};
            % Put frequently used parameters at the top level
            if ismember(currField, topFields)
                mD(iTrial).(currField) = expMetadata(iTrial).(currField);            
            else
                % Put the rest down a level in a separate struct
                mD(iTrial).expMetadata.(currField) = expMetadata(iTrial).(currField);
            end            
        end
        mD(iTrial).displayRate = double(mD(iTrial).displayRate);
    end
      
    % Extracted size and timing values
    mD(iTrial).volumeRate =  mD(iTrial).SI.hRoiManager.scanVolumeRate;
    mD(iTrial).nVolumes = mD(iTrial).SI.hFastZ.numVolumes;
    mD(iTrial).nPlanes = size(mD(iTrial).refImages, 3);
    mD(iTrial).volTimes = (1:1:mD(iTrial).nVolumes) / mD(iTrial).volumeRate;
    mD(iTrial).nDaqSamples = size(mD(iTrial).panelsPosData, 1);
    daqSampDur = mD(iTrial).nDaqSamples / mD(iTrial).trialDuration;
    mD(iTrial).daqSampTimes = (1:1:mD(iTrial).nDaqSamples) / mD(iTrial).SAMPLING_RATE;
    
    % Panels information
    xPosFunc = mD(iTrial).expMetadata.xDimPosFun.func;
    yPosFunc = mD(iTrial).expMetadata.yDimPosFun.func;
    mD(iTrial).panelsCycleFrames = numel((mD(iTrial).expMetadata.xDimPosFun.func));
    mD(iTrial).panelsCycleTime = mD(iTrial).panelsCycleFrames / mD(iTrial).displayRate;
    xPosRep = repmat(xPosFunc, 1, ceil(mD(iTrial).trialDuration / mD(iTrial).panelsCycleTime));
    yPosRep = repmat(yPosFunc, 1, ceil(mD(iTrial).trialDuration / mD(iTrial).panelsCycleTime));
    mD(iTrial).nPanelsFrames = mD(iTrial).displayRate * trialDuration;
    mD(iTrial).panelsPosX = xPosRep(1:mD(iTrial).nPanelsFrames);
    mD(iTrial).panelsPosY = yPosRep(1:mD(iTrial).nPanelsFrames);
    mD(iTrial).panelsFrameTimes = (1:1:mD(iTrial).nPanelsFrames) / mD(iTrial).displayRate;
    mD(iTrial).panelsPattern = mD(iTrial).expMetadata.pattern.Pats; % --> [y, x, frame, dim]
    
    % Opto stim timing
    if mD(iTrial).usingOptoStim
        
        % Make a list of each stim onset and offset throughout the entire trial
        stimTimes = mD(iTrial).optoStimTiming;
        mD(iTrial).optoStimOnsetTimes = (stimTimes(1):sum(stimTimes(2:3)):mD(iTrial).trialDuration ...
                - stimTimes(2))';
        mD(iTrial).optoStimOffsetTimes = (sum(stimTimes(1:2)):sum(stimTimes(2:3)) ...
            :mD(iTrial).trialDuration)';
    else
        mD(iTrial).optoStimOnsetTimes = [];
        mD(iTrial).optoStimOnsetTimes = [];
    end
    
end%iTrial
    
% Rename a couple of fields
[mD.daqSampRate] = mD.SAMPLING_RATE;
[mD.panelsDisplayRate] = mD.displayRate; 
mD = rmfield(mD, {'SAMPLING_RATE', 'displayRate'});

% Sort fields alphabetically
mD = orderfields(mD);

% Save data file
disp('Saving file...')
save(fullfile(parentDir, 'analysis_data.mat'), 'mD', '-v7.3');
disp('Saving complete');

