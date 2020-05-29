parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';

expList = load_expList('group', 'singleTrialAcq_with_ROIs')

% Load exp  and trial metadata
[expMd, trialMd] = load_metadata(expList);

% Load ROI data
roiData = load_roi_data(expList);

% Load event data
events = load_event_data(expList);

% Load fictrac data
ftData = load_ft_data(expList);

%% 

eventData = events.odor.eventData;
locEvents = events.locomotion;

% eventData = eventData(strcmp(eventData.concentration, 'neat'), :);
% eventData = eventData(strcmp(eventData.odorName, 'ACV'), :);

% Join ROI data, event data, and metadata tables
allEventData = inner_join(trialMd(:, {'expID', 'trialNum', 'nVolumes', 'originalTrialCount'}), ...
                            expMd(:, {'expID', 'volumeRate'}), roiData(:, {'expID', 'trialNum', ...
                            'roiName', 'rawFl', 'trialBaseline', 'expBaseline'}), [eventData, ...
                            table((1:size(eventData, 1))', 'VariableNames', {'eventNum'})]); 
                        
% Add FicTrac data, if applicable
if ~isempty(ftData)
    allEventData = innerjoin(allEventData, ftData(:, {'expID', 'trialNum', 'yawSpeed', ...
        'moveSpeed', 'frameTimes'}));
end

% Locomotion array
locEventArray = locEvents.create_logical_array(allEventData.nVolumes(1), ...
        allEventData.volumeRate(1));
locEventTrialList = locEvents.trialList;
locEventTrialList.locArrayRow = (1:size(locEventTrialList, 1))';
locEventTrialList = outerjoin(allEventData(:, {'expID', 'trialNum'}), locEventTrialList, 'type', 'left');
locEventTrialList.expID = locEventTrialList.expID_left;
locEventTrialList.trialNum = locEventTrialList.trialNum_left;
locEventTrialList = unique(locEventTrialList);

%%

% Set analysis window size in seconds
preEventWin = 4;
postEventWin = 4;

alignedFlData = [];
for iEvent = 1:size(eventData, 1)
    
    disp(['Processing event ', num2str(iEvent), ' of ', num2str(size(eventData, 1))]);
    
    currEventData = allEventData(allEventData.eventNum == iEvent, :);
    
    % Save the start and end times of the analysis window for the current event
    startTime = currEventData.onsetTime(1) - preEventWin;
    endTime = currEventData.onsetTime(1) + postEventWin;
    currEventData.startTime = ones(size(currEventData, 1), 1) * startTime;
    currEventData.endTime = ones(size(currEventData, 1), 1) * endTime;
    
    % Calculate sizes of the analysis window in volumes
    preEventWinVols = ceil(preEventWin * currEventData.volumeRate(1));
    postEventWinVols = ceil(postEventWin * currEventData.volumeRate(1));
    winSizeVols = 1 + preEventWinVols + postEventWinVols;
    
    % Calculate start and end volumes 
    onsetVol = round(currEventData.onsetTime(1) * currEventData.volumeRate(1));
    startVol = onsetVol - preEventWinVols;
    endVol = onsetVol + postEventWinVols;    
    
    % Extract raw fluorescence data in the event window
    eventFl = []; 
    eventLocVols = [];
    startPad = 0;
    endPad = 0;
    for iRoi = 1:size(currEventData, 1)        
        if startVol < 1
            startPad = -startVol + 1;
            startVol = 1;
        end
        if endVol > numel(currEventData.rawFl{iRoi})
            endPad = endVol - numel(currEventData.rawFl{iRoi});
            endVol = numel(currEventData.rawFl{iRoi});
        end
        eventFl{iRoi, 1} = [nan(startPad, 1); currEventData.rawFl{iRoi}(startVol:endVol); ...
            nan(endPad, 1)];
        
        % Extract locomotion volumes from logical array
        locRow = locEventTrialList(strcmp(locEventTrialList.expID, currEventData.expID{1}) & ...
                (locEventTrialList.trialNum == currEventData.trialNum(1)), {'locArrayRow'});
        if ~isnan(locRow.locArrayRow)
            currLocData = locEventArray(:, locRow.locArrayRow);
            eventLocVols{iRoi, 1} = [nan(startPad, 1); currLocData(startVol:endVol); ...
                    nan(endPad, 1)];
        else
            eventLocVols{iRoi, 1} = [nan(startPad, 1); 0 * (startVol:endVol)'; ...
                    nan(endPad, 1)]; 
        end
    end
    
    
    ftFrameTimes = currEventData.frameTimes{1};
    if ~isempty(ftFrameTimes)
        
        % Calculate average frame rate of FicTrac data
        avgFrameDur = ftFrameTimes(end) / numel(ftFrameTimes);
        avgFrameRate = 1 / avgFrameDur;
        
        % Resample FicTrac data at a constant frame rate
        targetSampleTimes = ftFrameTimes(1):avgFrameDur:ftFrameTimes(end);
        ftFrameInds = interp1(ftFrameTimes, 1:numel(ftFrameTimes), targetSampleTimes, 'nearest');
        ftFrameTimesResamp = ftFrameTimes(ftFrameInds);
        yawSpeedResamp = currEventData.yawSpeed{1}(ftFrameInds);
        moveSpeedResamp = currEventData.moveSpeed{1}(ftFrameInds);
        
        % Calculate sizes of the analysis window in resampled FicTrac frames
        preEventWinFrames = ceil(preEventWin * avgFrameRate);
        postEventWinFrames = ceil(postEventWin * avgFrameRate);
        winSizeFrames = 1 + preEventWinFrames + postEventWinFrames;
        
        % Calculate start and end frames
        onsetFrame = round(currEventData.onsetTime(1) * avgFrameRate);
        startFrame = onsetFrame - preEventWinFrames;
        endFrame = onsetFrame + postEventWinFrames;
        
        % Extract FicTrac data in the event window
        startPad = 0;
        endPad = 0;
        if startFrame < 1
            startPad = -startFrame + 1;
            startFrame = 1;
        end
        if endFrame > numel(ftFrameTimesResamp)
            endPad = endFrame - numel(ftFrameTimesResamp);
            endFrame = numel(ftFrameTimesResamp);
        end
        eventMoveSpeed = [nan(startPad, 1); moveSpeedResamp(startFrame:endFrame); ...
            nan(endPad, 1)];
        eventYawSpeed = [nan(startPad, 1); yawSpeedResamp(startFrame:endFrame); ...
            nan(endPad, 1)];
        eventFrameTimes = [nan(startPad, 1); ftFrameTimesResamp(startFrame:endFrame); ...
            nan(endPad, 1)];
        
        % Add data to currEvent table
        currEventData.eventFl = eventFl;
        currEventData.eventLocVols = eventLocVols;
        currEventData.eventMoveSpeed = repmat({eventMoveSpeed}, size(currEventData, 1), 1);
        currEventData.eventYawSpeed = repmat({eventYawSpeed}, size(currEventData, 1), 1);
        currEventData.eventFrameTimes = repmat({eventFrameTimes}, size(currEventData, 1), 1);
        
        % Remove original data fields and append to aligned data table
        alignedFlData = [alignedFlData; removevars(currEventData, {'rawFl', 'yawSpeed', 'moveSpeed', ...
            'frameTimes'})];
    
    end
end

%%

% test = alignedFlData(strcmp(alignedFlData.roiName, 'TypeD-L'), :);
test = alignedFlData(contains(alignedFlData.roiName, 'VLP-AMMC'), :);

flData = cell2mat(test.eventFl');
moveSpeed = cell2mat(test.eventMoveSpeed');
moveVols = cell2mat(test.eventLocVols');

trialBaseline = repmat(test.trialBaseline', size(flData, 1), 1);
expBaseline = repmat(test.expBaseline', size(flData, 1), 1);

dffData = (flData - trialBaseline) ./ trialBaseline;
expDffData = (flData - expBaseline) ./ expBaseline;

plotData = flData;
plotData = dffData;
plotData = expDffData;

nanTrials = sum(isnan(plotData), 1);
plotData = plotData(:, ~nanTrials);
moveVols = moveVols(:, ~nanTrials);
moveSpeed = smoothdata(moveSpeed(:, ~nanTrials), 1, 'gaussian', 5);
moveSpeedSorted = sort(moveSpeed(:));
capVal = moveSpeedSorted(round(numel(moveSpeedSorted) * 0.999));
moveSpeed(moveSpeed > capVal) = capVal;

moveVols(1, 134) = 1;

plotData = plotData(:, ~any(moveVols, 1));
moveSpeed = moveSpeed(:, ~any(moveVols, 1));
moveVols = moveVols(:, ~any(moveVols, 1));



figure(1);clf; imagesc(plotData');
title('Fl')

figure(3);clf; imagesc(moveSpeed');
title('Move Speed')

figure(4);clf; imagesc(moveVols');
title('Locomotion')

figure(2);clf; 
% plot(smoothdata(plotData, 1, 'gaussian', 5))
hold on; plot(mean(plotData, 2, 'omitnan'), 'linewidth', 3, 'color', 'k')

% Shade +/- SEM
sd = std(plotData, [], 2, 'omitnan');
avg = mean(smoothdata(plotData, 'gaussian', 3, 'omitnan'), 2);
sem = sd ./ (size(plotData, 2)^0.5);
jbfill(1:size(plotData, 1), [avg + sem]', [avg - sem]', 'k', 'k', 1, 0.2);

disp(unique(alignedFlData.roiName))

