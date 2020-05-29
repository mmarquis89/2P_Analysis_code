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
                if ~isempty(shutoffVols)
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
                dataTable = generate_data_table(obj, loadOneNoteData, oneNoteDataFile);
            else
                dataTable = generate_data_table(obj, loadOneNoteData);
            end
            
            % Get list of all field names
            allFieldNames = dataTable.Properties.VariableNames;
            
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
        
        % OUTPUT A TABLE WITH PART OR ALL OF THE DATA
        function dataTable = output_analysis_subset(obj, analysisWin, filterDefs)
            
            % Generate table
            dataTable = generate_data_table(obj);
            
            % Split the time-vector-based fields off of filterDefs so that the rest of the filters
            % can be applied before they are trimmed to the size of the analysis window
            if nargin > 2
                eventNames = fieldnames(obj.eventObjects);
                firstRoundFilterDefs = rmfield(filterDefs, [eventNames; {'moveSpeed'; 'yawSpeed'}]);
            end
            
            % Subset table if filterDefs were provided 
            if nargin > 2
                dataTable = subset_data_table(dataTable, firstRoundFilterDefs);
            end
                
            if ~isempty(dataTable)
                
                % Trim to analysis window if one was provided
                if nargin > 1 && ~isempty(analysisWin)
                    dataTable = trim_data_table(dataTable, analysisWin, fieldnames(obj.eventObjects));
                end
                
                % Second round of subsetting for time-vector-based fields
                if nargin > 2
                    dataTable = subset_data_table(dataTable, filterDefs);
                end
            end
            
            % Save a copy of the filterDefs and analysisWin to the table's custom properties
            dataTable = addprop(dataTable, {'analysisWin', 'filterDefs'}, {'table', 'table'});
            dataTable.Properties.CustomProperties.analysisWin = analysisWin;
            dataTable.Properties.CustomProperties.filterDefs = filterDefs;
        end

    end%Methods (ordinary)
    
    methods(Static)
        
    end%Methods (static)
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
function dataTable = generate_data_table(obj, loadOneNoteData, oneNoteDataFile)

if nargin < 2
    loadOneNoteData = 1;
end
if nargin < 3 && loadOneNoteData
    oneNoteDataFile = ['D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\', ...
            'OneNoteMetadata.csv'];
end

% Create base table
dataTable = inner_join(obj.sourceMd, obj.alignEventTable, obj.eventFlTable, ...
    obj.eventFilterTable);

% Load and join OneNote metadata if necessary
if loadOneNoteData
    oneNoteData = readtable(oneNoteDataFile, 'delimiter', ',');
    dataTable = inner_join(dataTable, oneNoteData);
end

% Add a custom property to keep track of the alignment event name
dataTable = addprop(dataTable, 'alignEventName', 'table');
dataTable.Properties.CustomProperties.alignEventName = obj.alignEventName;
            
end

% Trim all aligned data fields in a data table to a specified time window
function trimmedTable = trim_data_table(dataTable, analysisWin, eventNames)

volTimes = {}; eventFl = {}; eventData = struct();
frameTimes = {}; moveSpeed = {}; yawSpeed = {};
for iRow = 1:size(dataTable, 1)
    if ~mod(iRow, 1000)
        disp([num2str(iRow), ' of ', num2str(size(dataTable, 1))]);
    end
    currRow = dataTable(iRow, :);
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
trimmedTable = dataTable;
trimmedTable.volTimes = volTimes;
trimmedTable.eventFl = eventFl;
for iType = 1:numel(eventNames)
    trimmedTable.(eventNames{iType}) = {eventData.(eventNames{iType})}';
end
trimmedTable.frameTimes = frameTimes;
trimmedTable.moveSpeed = moveSpeed;
trimmedTable.yawSpeed = yawSpeed;

end%function

% Use filterDefs to return a subset of rows from the full data table
function outputTable = subset_data_table(dataTable, filterDefs)
%  If filter is a string:
%   - Apply a regex string to a character field
%       Example: 'regex'
% 
%   - Evaluate as an expression to apply to a numeric field after replacing 'x' with 'filtData'
%       Example: '(x > 4 | x == 2) & x ~= 20'
% 
% If filter is numeric:
%   - Select a specific value or set of values from a numeric scalar field:
%       Example: 2 or [2 4 6]
% 
%   - Specify a code for a logical condition to apply one of the logical event vector fields:
%       0:  remove events with any occurance of the filter event
%       -1: remove events with any pre-onset occurance of the filter event
%       -2: remove events with any post-onset occurance of the filter event
%       1: remove events withOUT any pre-onset occurance of the filter event
%       2: remove events withOUT any post-onset occurance of the filter event

filterVec = ones(size(dataTable, 1), 1);
filtNames = fieldnames(filterDefs);
alignEventName = dataTable.Properties.CustomProperties.alignEventName;
for iFilt = 1:numel(filtNames)
    
    % Get current filter info
    filtName = filtNames{iFilt};
    filtValue = filterDefs.(filtName);
    filtData = dataTable.(filtName);
    
    % Determine the data type of the filter value and the table field to be filtered
    filterType = get_filter_type(filtValue);
    dataType = get_data_type(filtData);
    
    % ------ Decide how to interpret the filter ------
    if ~isempty(filtValue)
        
        % A regex string for a char field
        if strcmp(filterType, 'charVector') && strcmp(dataType, 'charVector')            
            currFiltVec = ~cellfun(@isempty, regexp(filtData, filtValue, 'once'));
            
        % A string to be evaluated as an expression after replacing 'x' with 'filtData' 
        elseif strcmp(filterType, 'charVector')
            currFiltVec = eval(regexprep(filtValue, 'x', 'filtData'));
        
        % A value specifying a logical filter vector condition
        elseif strcmp(filterType, 'numericScalar') && strcmp(dataType, 'numericVector')
            
             onsetVols = cellfun(@(x) argmin(abs(x)), dataTable.volTimes);
             onsetVols = mat2cell(onsetVols, ones(size(onsetVols)));
             
            % Filter out events with a pre-onset occurance of the alignment event type
            if filtValue < 1 && strcmp(filtName, alignEventName)
               
                currFiltVec = ~cellfun(@(x, y) any(x(1:(y - 1))), filtData, onsetVols);
            
            % Filter out events withOUT a pre-onset occurance of the alignment event time
            elseif filtValue == 1 && strcmp(filtName, alignEventName)

                currFiltVec = cellfun(@(x, y) any(x(1:(y - 1))), filtData, onsetVols);                
                
            % Filter out events with any occurance of the filter event type
            elseif filtValue == 0
                currFiltVec = ~cellfun(@any, filtData);
            
            % Filter out events with any pre-onset occurance of the filter event type
            elseif filtValue == -1
                currFiltVec = ~cellfun(@(x, y) any(x(1:(y - 1))), filtData, onsetVols);
                
            % Filter out events withOUT any pre-onset occurance of the filter event type
            elseif filtValue == 1
                currFiltVec = cellfun(@(x, y) any(x(1:(y - 1))), filtData, onsetVols);
                
            % Filter out events with any post-onset occurance of the filter event type
            elseif filtValue == -2
                currFiltVec = ~cellfun(@(x, y) any(x(y:end)), filtData, onsetVols);
                
            % Filter out events withOUT any post-onset occurance of the filter event type
            elseif filtValue == 2
                currFiltVec = cellfun(@(x, y) any(x((y:end))), filtData, onsetVols);                
                        
            else
                error(['Invalid filter type "', filterType, '" for field ', filtName]); 
            end
            
        % A specific scalar value or set of acceptable scalar values
        elseif strcmp(dataType, 'numericScalar')            
            currFiltVec = ismember(filtData, filtValue);    
            
        else
            error(['Invalid filter type "', filterType, '" for field ', filtName]); 
        end
        
        filterVec = filterVec .* currFiltVec;
    end    
end%iFilt

% Subset data table
filterVec = logical(filterVec);
outputTable = dataTable(filterVec, :);

end%function

% Determine data type of a particular filterDef field value
function filterDataType = get_filter_type(filtValue)
if isa(filtValue, 'char')
    filterDataType = 'charVector';
elseif isa(filtValue, 'numeric') && isscalar(filtValue)
    filterDataType = 'numericScalar';
else
    filterDataType = 'numericVector';
end
end

% Determine data type of a particular table field to use in filtering
function fieldDataType = get_data_type(filtData)
if isa(filtData, 'cell') && strcmp(unique(cellfun(@class, filtData, 'uniformOutput', 0)), ...
        'char')
    fieldDataType = 'charVector';
elseif isa(filtData, 'cell') && strcmp(unique(cellfun(@class, filtData, 'uniformOutput', 0)), ...
        'double')
    fieldDataType = 'numericVector';
else
    fieldDataType = 'numericScalar';
end
end


