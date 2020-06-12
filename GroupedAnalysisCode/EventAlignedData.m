classdef EventAlignedData
% ==================================================================================================   
%    
% Properties:
%       dataDir
%       alignEventName
%       maxEventWin
%       eventObjects
%       sourceMd
%       alignEventTable
%       eventFlTable 
%       eventFilterTable
%
% Methods:
%       obj = set_align_event(eventName)
%       obj = extract_roi_data()
%       obj = extract_fictrac_data()
%       obj = create_filter_event_vectors()
%       filterDefs = create_filterDefs(loadOneNoteData, oneNoteDataFile, removeFieldsOverride)
%       dataTable = output_analysis_subset(analysisWin, filterDefs)
%
% Subclasses:
%
%
% ==================================================================================================
    properties
        dataDir
        
        alignEventName
        maxEventWin
        eventObjects % Struct with a field containing an event object for each event time        
        sourceMd % table with originally loaded data
        
        alignEventTable
        eventFlTable
        eventFilterTable
    end
    
    properties (Dependent)
        trialList
    end
    
    methods
        
        % CONSTRUCTOR
        function obj = EventAlignedData(expList, varargin) 
            
            p = inputParser;
            addParameter(p, 'dataDir', ...
                    'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments');
            addParameter(p, 'loadFicTrac', 1);
            parse(p, varargin{:});
            dataDir = p.Results.dataDir;
            obj.dataDir = dataDir;
            
            % Load all exp and trial metadata files
            [expMd, trialMd] = load_metadata(expList, dataDir);
            
            % Load ALL event data objects
            obj.eventObjects = load_event_data(expList, dataDir);
            
            % Create source data table
            obj.sourceMd = inner_join(trialMd(:, {'expID', 'trialNum', 'nVolumes', 'trialDuration', ...
                        'originalTrialCount', 'pmtShutoffVols'}), expMd(:, {'expID', 'expName', ...
                        'volumeRate'}));
        end        
        
        % GET/SET METHODS FOR DEPENDENT PROPERTIES
        function trialList = get.trialList(obj)
            trialList = unique(obj.sourceMd(:, {'expID', 'trialNum'}));
        end
        
        % SET ALIGNMENT EVENT
        function obj = set_align_event(obj, eventName)
            
            obj.maxEventWin = [10, 10];
            
            obj.alignEventName = eventName;
            alignEventObj = obj.eventObjects.(eventName);
            
            % Add some columns from the sourceData table, only using unique rows when joining ...
            % to the alignEventData
            alignEventData = innerjoin(alignEventObj.eventData, unique(obj.sourceMd(:, {'expID', ...
                    'trialNum', 'nVolumes', 'volumeRate', 'trialDuration', 'originalTrialCount'})));
            
            % Calculate indices of aligned sample data and necessary nan padding size
            disp('Aligning data to primary event type...')
            tic
            startTime = []; endTime = [];
            startPadVols = []; endPadVols = [];
            onsetVol = []; startVol = []; endVol = [];
            for iEvent = 1:size(alignEventData, 1)   
                if ~mod(iEvent, 1000)
                    disp([num2str(iEvent), ' of ', num2str(size(alignEventData, 1))]);
                end
                currEventData = alignEventData(iEvent, :);
                startTime(iEvent) = currEventData.onsetTime(1) - obj.maxEventWin(1);
                endTime(iEvent) = currEventData.onsetTime(1) + obj.maxEventWin(2);
                
                % Calculate sizes of the analysis window in samples
                preEventWinVols = ceil(obj.maxEventWin(1) * currEventData.volumeRate);
                postEventWinVols = ceil(obj.maxEventWin(2) * currEventData.volumeRate);
                
                % Calculate onset, start, and end volumes
                expDateNum = str2double(currEventData.expID{1}(1:8));
%                 disp(expDateNum)
                if expDateNum > 20181010 && expDateNum < 20200101 
                    % Alternate volTimes calculation for Gapless acq experiments
                    volsPerTrial = currEventData.nVolumes / currEventData.originalTrialCount;
                    trialVolTimes = (1:volsPerTrial) / currEventData.volumeRate;
                    volTimes = [];
                    for iTrial = 1:currEventData.originalTrialCount
                        volTimes = [volTimes, (trialVolTimes + ((currEventData.trialDuration / ...
                                currEventData.originalTrialCount) * (iTrial - 1)))];
                    end
                    onsetVol(iEvent) = argmin(abs(volTimes - currEventData.onsetTime));
                else
                    onsetVol(iEvent) = round(currEventData.onsetTime * currEventData.volumeRate);
                end
                startVol(iEvent) = onsetVol(iEvent) - preEventWinVols;
                endVol(iEvent) = onsetVol(iEvent) + postEventWinVols;
                
                % Calculate size of nan padding
                if startVol(iEvent) < 1
                    startPadVols(iEvent) = -startVol(iEvent) + 1;
                    startVol(iEvent) = 1;
                else
                    startPadVols(iEvent) = 0;
                end
                if endVol(iEvent) > currEventData.nVolumes
                    endPadVols(iEvent) = endVol(iEvent)- currEventData.nVolumes;
                    endVol(iEvent) = currEventData.nVolumes;
                else
                    endPadVols(iEvent) = 0;
                end
                                    
            end
            
            % Create alignEventTable
            obj.alignEventTable = [removevars(alignEventData, {'nVolumes', 'volumeRate'}), ...
                    table(startTime', endTime', startPadVols', endPadVols', onsetVol', startVol', endVol', ...
                    'VariableNames', {'startTime', 'endTime', 'startPadVols', 'endPadVols', 'onsetVol', ...
                    'startVol', 'endVol'})];
            disp(['Alignment completed in ', num2str(round(toc)), ' seconds']);
        end
        
        % EXTRACT ROI DATA
        function obj = extract_roi_data(obj)
            
            % Load ROI data
            expList = unique(obj.trialList(:, 1));
            roiData = load_roi_data(expList, obj.dataDir);
            
            % Create table of imaging data
            sourceFlData = inner_join(obj.sourceMd(:, {'expID', 'trialNum', 'pmtShutoffVols'}), ...
                    roiData(:, {'expID', 'trialNum', 'roiName', 'rawFl', 'trialBaseline', ...
                    'expBaseline'}));
                
            % Set any volumes when the PMT was turned off to 'nan'
            for iTrial = 1:size(sourceFlData, 1)
                shutoffVols = sourceFlData.pmtShutoffVols{iTrial};
                if ~isempty(shutoffVols) & ~all(isnan(shutoffVols))
                    sourceFlData.rawFl{iTrial}(shutoffVols) = nan;
                end
            end
            
            % Join with alignment event table 
            obj.eventFlTable = inner_join(obj.sourceMd(:, {'expID', 'trialNum'}), ...
                    sourceFlData, obj.alignEventTable);
                
            % Extract fluorecence data in analysis window around each event
            disp('Extracting fluorescence data for each event...')
            eventFl = {};
            for iRow = 1:size(obj.eventFlTable, 1)
                if ~mod(iRow, 1000)
                    disp([num2str(iRow), ' of ', num2str(size(obj.eventFlTable, 1))]);
                end
                currRow = obj.eventFlTable(iRow, :);
                eventFl{iRow, 1} =  [nan(currRow.startPadVols, 1); currRow.rawFl{1}( ...
                        currRow.startVol:currRow.endVol); nan(currRow.endPadVols, 1)];
            end
            obj.eventFlTable.eventFl = eventFl;
            obj.eventFlTable = obj.eventFlTable(:, {'expID', 'trialNum', 'onsetTime', 'roiName', ...
                    'eventFl', 'trialBaseline', 'expBaseline'});
            disp('Fluorescence data extracted.')
        end
        
        % EXTRACT FICTRAC DATA
        function obj = extract_fictrac_data(obj)
            
            % Load FicTrac data if applicable
            expList = unique(obj.trialList(:, 1));
            ftData = load_ft_data(expList, obj.dataDir);
            
            % Add ficTracData if there is any
            if ~isempty(ftData)
                
                % Remove any existing FicTrac fields from the sourceMd table
                obj.sourceMd = remove_fictrac_fields(obj.sourceMd);
                obj.alignEventTable = remove_fictrac_fields(obj.alignEventTable);
                
                % Add newly-loaded FicTrac data to table
                ftData = outerjoin(obj.sourceMd, ftData(:, {'expID', 'trialNum', ...
                        'moveSpeed', 'yawSpeed', 'frameTimes'}), 'type', 'left', 'mergekeys', 1);

                eventFt = inner_join(ftData, obj.alignEventTable);
                disp('Extracting FicTrac data within event windows...')
                emptyVec = nan(1, size(eventFt, 1));
                emptyCell = repmat({nan}, size(eventFt, 1), 1);
                startPadFrames = emptyVec; endPadFrames = emptyVec;
                onsetFrame = emptyVec; startFrame = emptyVec; endFrame = emptyVec; 
                avgFrameRate = emptyVec;
                ftFrameTimesResamp = emptyCell; 
                moveSpeedResamp = emptyCell; 
                yawSpeedResamp = emptyCell;
                for iRow = 1:size(eventFt, 1)
                    if ~mod(iRow, 1000)
                        disp([num2str(iRow), ' of ', num2str(size(eventFt, 1))]);
                    end
                    currRow = eventFt(iRow, :);
                    ftFrameTimes = currRow.frameTimes{1};
                    if ~isempty(ftFrameTimes)

                        % Calculate average frame rate of FicTrac data
                        avgFrameDur = ftFrameTimes(end) / numel(ftFrameTimes);
                        avgFrameRate(iRow) = 1 / avgFrameDur;

                        % Resample FicTrac data at a constant frame rate
                        targetSampleTimes = ftFrameTimes(1):avgFrameDur:ftFrameTimes(end);
                        ftFrameInds = interp1(ftFrameTimes, 1:numel(ftFrameTimes), targetSampleTimes, ...
                            'nearest');
                        ftFrameTimesResamp{iRow} = ftFrameTimes(ftFrameInds);
                        moveSpeedResamp{iRow} = currRow.moveSpeed{1}(ftFrameInds);
                        yawSpeedResamp{iRow} = currRow.yawSpeed{1}(ftFrameInds);

                        %Calculate sizes of the analysis window in frames
                        preEventWinFrames = ceil(obj.maxEventWin(1) * avgFrameRate(iRow));
                        postEventWinFrames = ceil(obj.maxEventWin(2) * avgFrameRate(iRow));
                        
                        % Calculate onset, start, and end frames
                        onsetFrame(iRow) = round(currRow.onsetTime * avgFrameRate(iRow));
                        startFrame(iRow) = onsetFrame(iRow) - preEventWinFrames;
                        endFrame(iRow) = onsetFrame(iRow) + postEventWinFrames;                         
                        
                        % Calculate size of nan padding
                        if startFrame(iRow) < 1
                            startPadFrames(iRow) = -startFrame(iRow) + 1;
                            startFrame(iRow) = 1;
                        else
                            startPadFrames(iRow) = 0;
                        end
                        if endFrame(iRow) > numel(ftFrameTimesResamp{iRow})
                            endPadFrames(iRow) = endFrame(iRow) - numel(ftFrameTimesResamp{iRow});
                            endFrame(iRow) = numel(ftFrameTimesResamp{iRow});
                        else
                            endPadFrames(iRow) = 0;
                        end
                    end%if
                end%iRow
                eventMoveSpeed = emptyCell; eventYawSpeed = emptyCell; eventFrameTimes = emptyCell;
                for iRow = 1:size(eventFt, 1)
                    if ~isempty(eventFt.frameTimes{iRow})
                        eventMoveSpeed{iRow, 1} = [nan(startPadFrames(iRow), 1); moveSpeedResamp{iRow}( ...
                            startFrame(iRow):endFrame(iRow)); nan(endPadFrames(iRow), 1)];
                        eventYawSpeed{iRow, 1} = [nan(startPadFrames(iRow), 1); yawSpeedResamp{iRow}( ...
                            startFrame(iRow):endFrame(iRow)); nan(endPadFrames(iRow), 1)];
                        eventFrameTimes{iRow, 1} = [nan(startPadFrames(iRow), 1); ftFrameTimesResamp{iRow}( ...
                            startFrame(iRow):endFrame(iRow)); nan(endPadFrames(iRow), 1)];
                    end
                end%iRow

                eventFtTable = eventFt(:, {'expID', 'trialNum', 'onsetTime'});
                eventFtTable.startPadFrames = startPadFrames';
                eventFtTable.endPadFrames = endPadFrames';
                eventFtTable.onsetFrame = onsetFrame';
                eventFtTable.startFrame = startFrame';
                eventFtTable.endFrame = endFrame';
                eventFtTable.moveSpeed = eventMoveSpeed;
                eventFtTable.yawSpeed = eventYawSpeed;
                eventFtTable.frameTimes = eventFrameTimes;
                eventFtTable.avgFrameRate = avgFrameRate';

                % Remove any existing FicTrac fields from the alignEvent table
                obj.alignEventTable = remove_fictrac_fields(obj.alignEventTable);
                
                obj.alignEventTable = outerjoin(obj.alignEventTable, eventFtTable, 'type', 'left', ...
                    'mergekeys', 1);
                disp('FicTrac extraction complete.')
            end%if
        end%function
        
        % CREATE FILTER EVENT VECTORS
        function obj = create_filter_event_vectors(obj)
            filterTable = [];
            eventObjNames = fieldnames(obj.eventObjects);
            disp('Creating filter event vectors...')
            tic
            
            % Create logical arrays for each event type
            disp('Creating logical arrays...')
            for iObj = 1:numel(eventObjNames)
                
                currObjName = eventObjNames{iObj};
                currEventObj = obj.eventObjects.(currObjName);
                filterEventData = innerjoin(currEventObj.eventData, obj.sourceMd(:, {'expID', ...
                        'trialNum', 'nVolumes', 'volumeRate', 'originalTrialCount', ...
                        'trialDuration'}));
                disp(currObjName)
                trialList = unique(filterEventData(:, {'expID', 'trialNum', 'nVolumes', ...
                        'originalTrialCount', 'trialDuration', 'volumeRate'}));
                    
                % Don't bother processing trials in which the alignment event never occurs
                trialList = innerjoin(trialList, unique(obj.alignEventTable(:, ...
                        {'expID', 'trialNum'})));
                    
                for iTrial = 1:size(trialList, 1)
                    if ~mod(iTrial, 1000)
                        disp([num2str(iTrial), ' of ', num2str(size(trialList, 1))]);
                    end
                    currTrial = trialList(iTrial, :);
                    
                    % Calculate volume times
                    volTimes = calc_volTimes(currTrial.nVolumes, currTrial.volumeRate, ...
                            currTrial.trialDuration, currTrial.originalTrialCount);
                    
                    % Get logical array for current event and trial volumes
                    eventArr = currEventObj.create_logical_array(currTrial.nVolumes, ...
                            volTimes, currTrial);
                    newRow = [currTrial(:, {'expID', 'trialNum'}), table(mat2cell(eventArr, ...
                            size(eventArr, 1), ones(size(eventArr, 2), 1)), repmat({currObjName}, ...
                            size(eventArr, 2), 1), 'variableNames', {'eventVector', 'eventName'})];
                    filterTable = [filterTable; newRow];
                end
            end
            filterVectorTable = filterTable;
            
            % Extract aligned data from logical arrays
            disp('Aligning data from logical arrays...')
            obj.eventFilterTable = obj.alignEventTable;
            for iType = 1:numel(eventObjNames)
                currName = eventObjNames{iType};
                disp(currName)
                currEventData = outerjoin(obj.alignEventTable, filterVectorTable(...
                        strcmp(filterVectorTable.eventName, currName), :), 'type', 'left', ...
                        'mergekeys', 1);
                currEventAlignedVectors = {};
                for iEvent = 1:size(currEventData, 1)
                    currRow = currEventData(iEvent, :);
                    if ~isempty(currRow.eventVector{:})
                        currEventAlignedVectors{iEvent, 1} = [nan(currRow.startPadVols, 1); ...
                            currRow.eventVector{:}(currRow.startVol:currRow.endVol); ...
                            nan(currRow.endPadVols, 1)];
                    else
                        currEventAlignedVectors{iEvent, 1} = [nan(currRow.startPadVols, 1); ...
                            zeros(numel(currRow.startVol:currRow.endVol), 1); ...
                            nan(currRow.endPadVols, 1)];
                    end
                end
                obj.eventFilterTable.(currName) = currEventAlignedVectors;
            end
            obj.eventFilterTable = obj.eventFilterTable(:, [{'expID', 'trialNum', 'onsetTime'}, ...
                    eventObjNames']);
            disp(['Filter event vectors created in ', num2str(round(toc)), ' seconds'])
        end
        
        % CREATE A FILTER DEFS OBJECT
        function filterDefs = create_filterDefs(obj, varargin)
           
            % Parse optional arguments
            p = inputParser;
            addParameter(p, 'loadOneNoteData', 1);
            addParameter(p, 'oneNoteDataFile', []);
            addParameter(p, 'removeFieldsOverride', 0);
            parse(p, varargin{:});
            loadOneNoteData = p.Results.loadOneNoteData;
            oneNoteDataFile = p.Results.oneNoteDataFile;
            removeFieldsOverride = p.Results.removeFieldsOverride;
            
            % Create a full data table
            if ~isempty(oneNoteDataFile)
                dataTable = generate_DataTable(obj, loadOneNoteData, oneNoteDataFile);
            else
                dataTable = generate_DataTable(obj, loadOneNoteData);
            end
            
            % Get list of all field names
            allFieldNames = dataTable.sourceData.Properties.VariableNames;
            
            % Remove names of fields that don't make sense to use as filters (doing this 
            % subtractively instead of selecting desired fields because table may have stimEvent
            % metadata fields with unknown names) 
            REMOVE_LIST = {'pmtShutoffVols', 'volumeRate', 'onsetTime', 'offsetTime', 'startTime', ...
                    'endTime', 'startPadVols', 'endPadVols' 'onsetVol', 'startVol', 'endVol', ...
                    'startPadFrames', 'endPadFrames', 'onsetFrame', 'startFrame', 'endFrame', ...
                    'frameTimes', 'avgFrameRate', 'eventFl', 'trialBaseline' 'expBaseline', ...
                    'ballFlowRate', 'purpose', 'prepNotes', 'volTimes'};         
            if isa(removeFieldsOverride, 'cell')
                REMOVE_LIST = REMOVE_LIST(~ismember(REMOVE_LIST, removeFieldsOverride));
            elseif removeFieldsOverride == 1
                REMOVE_LIST = {''};
            end
            filterFieldNames = allFieldNames(~ismember(allFieldNames, REMOVE_LIST));            
            
            % Create structure to hold filter definitions
            filterDefs = struct();
            for iField = 1:numel(filterFieldNames)
                filterDefs.(filterFieldNames{iField}) = [];
            end
            
            
        end
        
        % OUTPUT A DataTable OBJECT WITH PART OR ALL OF THE DATA
        function outputDataTable = output_analysis_DataTable(obj, analysisWin)
            
            % Generate table
            outputDataTable = generate_DataTable(obj);
            
            % Save a the analysisWin and alignEventName to the DataTable's custom properties
            outputDataTable.customProps.analysisWin = analysisWin;
            outputDataTable.customProps.alignEventName = obj.alignEventName;
                
            % Trim the DataTable's sourceData to the analysis window if one was provided
            if nargin > 1 && ~isempty(analysisWin)
                trimmedSourceData = trim_aligned_data(outputDataTable.sourceData, analysisWin, ...
                        fieldnames(obj.eventObjects));
                outputDataTable.sourceData = trimmedSourceData;
                outputDataTable.filterMat = true(size(trimmedSourceData));
            end
        end

    end%Methods (ordinary)
    
end%Classdef

% ==================================================================================================
% Local functions
% ==================================================================================================

% Clear any FicTrac-related fields from the input table to prevent errors when re-loading data
function outputTable = remove_fictrac_fields(inputTable)
ftFieldNames = {'moveSpeed', 'yawSpeed', 'frameTimes', 'startPadFrames', 'endPadFrames', ...
        'onsetFrame', 'startFrame', 'endFrame', 'eventMoveSpeed', 'eventYawSpeed', ...
        'eventFrameTimes'};

tableVarNames = inputTable.Properties.VariableNames;
existingFieldNames = ftFieldNames(contains(ftFieldNames, tableVarNames));
outputTable = removevars(inputTable, existingFieldNames);

end

% Generate comprehensive data table from all aligned data
function outputTable = generate_DataTable(obj, loadOneNoteData, oneNoteDataFile)

if nargin < 2
    loadOneNoteData = 1;
end
if nargin < 3 && loadOneNoteData
    oneNoteDataFile = ['D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\', ...
            'OneNoteMetadata.csv'];
end

% Create base table
fullTable = inner_join(obj.sourceMd, obj.alignEventTable, obj.eventFlTable, ...
    obj.eventFilterTable);

% Load and join OneNote metadata if necessary
if loadOneNoteData
    oneNoteData = readtable(oneNoteDataFile, 'delimiter', ',');
    fullTable = inner_join(fullTable, oneNoteData);
end

outputTable = DataTable(fullTable);

end

% Trim all aligned data fields in a table to a specified time window
function trimmedTable = trim_aligned_data(inputTable, analysisWin, eventNames)

volTimes = {}; eventFl = {}; eventData = struct();
frameTimes = {}; moveSpeed = {}; yawSpeed = {};
for iRow = 1:size(inputTable, 1)
    if ~mod(iRow, 1000)
        disp([num2str(iRow), ' of ', num2str(size(inputTable, 1))]);
    end
    currRow = inputTable(iRow, :);
    relOnsetTime = currRow.onsetTime - currRow.startTime;
    
    % Trim volume fields (fluorescence data and logical filter vectors)
    relVolTimes = ((1:numel(currRow.eventFl{1})) / currRow.volumeRate) ...
        - relOnsetTime;
    eventVolInds = relVolTimes > -analysisWin(1) & relVolTimes < analysisWin(2);
    
    volTimes{iRow, 1} = relVolTimes(eventVolInds);
    eventFl{iRow, 1} = currRow.eventFl{:}(eventVolInds);
    for iType = 1:numel(eventNames)
        eventData(iRow).(eventNames{iType}) = ...
            currRow.(eventNames{iType}){1}(eventVolInds);
    end
    
    % Trim FicTrac frame fields
    if ~isnan(currRow.avgFrameRate)
        relFrameTimes = ((1:numel(currRow.moveSpeed{1})) / currRow.avgFrameRate) ...
            - relOnsetTime;
        eventFrameInds = relFrameTimes > - analysisWin(1) & ...
            relFrameTimes < analysisWin(2);
        frameTimes{iRow, 1} = relFrameTimes(eventFrameInds)';
        moveSpeed{iRow, 1} = currRow.moveSpeed{:}(eventFrameInds);
        yawSpeed{iRow, 1} = currRow.yawSpeed{:}(eventFrameInds);
    else
        frameTimes{iRow, 1} = nan;
        moveSpeed{iRow, 1} = nan;
        yawSpeed{iRow, 1} = nan;
    end
end%iRow

% Overwrite original fields with trimmed data
trimmedTable = inputTable;
trimmedTable.volTimes = volTimes;
trimmedTable.eventFl = eventFl;
for iType = 1:numel(eventNames)
    trimmedTable.(eventNames{iType}) = {eventData.(eventNames{iType})}';
end
trimmedTable.frameTimes = frameTimes;
trimmedTable.moveSpeed = moveSpeed;
trimmedTable.yawSpeed = yawSpeed;

end%function





