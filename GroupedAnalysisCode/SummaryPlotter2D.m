classdef SummaryPlotter2D
% ==================================================================================================   
%    
% Properties:
%
%
%
% Methods:
%
%           
%
%
% ==================================================================================================

    properties
        parentDir
        sourceDataTable
        eventData
        params
        savedPlotParams
    end
    properties (SetAccess = immutable)
        FRAME_RATE
    end
    properties (Dependent)
        expList
    end
    
    methods
        
        % CONSTRUCTOR
        function obj = SummaryPlotter2D(expListInput)
            
            if nargin < 1
                % If called without any arguments, default to using full expList 
                expList = load_expList();
            elseif isa(expListInput, char) 
                % If input arg is a string, use it to filter the full expList by groupName
                expList = load_expList();
                expList = expList(contains(expList.expName, expListInput), :);
            else
                expList = expListInput;
            end
%             obj.expList = expList;
            obj.parentDir = ...
                    'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';
                
            obj.FRAME_RATE = 25;
            allExpData = [];
            for iExp = 1:size(expList, 1)

                currExpID = expList.expID{iExp};
                disp(currExpID);

                % Generate full paths to data files
                roiFileName = [currExpID, '_roiData.mat'];
                roiFile = fullfile(obj.parentDir, roiFileName);
                ftFileName = [currExpID, '_ficTracData.mat'];
                ftFile = fullfile(obj.parentDir, ftFileName);
                expMdFile = fullfile(obj.parentDir, [currExpID, '_expMetadata.csv']);
                trialMdFile = fullfile(obj.parentDir, [currExpID, '_trialMetadata.mat']);

                if exist(roiFile, 'file') && exist(expMdFile, 'file')
                    load(roiFile, 'roiData');
                    if exist(ftFile, 'file')
                        load(ftFile, 'ftData');
                    else
                        nanCol = repmat({nan}, size(roiData, 1), 1);
                        ftData = table(roiData.expID, roiData.trialNum, roiData.roiName, nanCol, ...
                                nanCol, nanCol, ...
                                nanCol, nanCol, nanCol, nanCol, nanCol, nanCol, nanCol, nanCol, ...
                                nanCol, 'variableNames', {'expID', 'trialNum', 'roiName', 'intX', ...
                                'intY', 'intHD', 'moveSpeed', 'intFwMove', 'intSideMove', 'yawSpeed', ...
                                'fwSpeed', 'sideSpeed', 'frameTimes', 'badVidFrames', 'meanFlow'});
                    end
                    expMd = readtable(expMdFile, 'delimiter', ',');
                    load(trialMdFile, 'trialMetadata');
                    trialMd = trialMetadata;
                    allExpData = [allExpData; inner_join(expMd, trialMd, roiData, ftData)];
                else
                    disp(['Skipping ', currExpID, ' due to one or more missing files']);
                end
            end%iExp
            
            % Create source DataTable
            obj.sourceDataTable = DataTable(inner_join(expList, allExpData));
            
            % Load event data
            obj.eventData = load_event_data(unique(obj.sourceDataTable.sourceData.expID));
            
            % Initialize default parameters
            obj.params = initialize_default_params(obj);
            obj.savedPlotParams = [];
        end
        
        % GET/SET METHODS FOR DEPENDENT PROPERTIES
        function expList = get.expList(obj)
            expList = obj.sourceDataTable.sourceData(:, {'expID', 'expName', 'groupName'});
        end
        
        % SAVE MAX AND MIN VALS FOR FUTURE USE
        function obj = save_colormap_range(obj)
            newRow = table({obj.params.expID}, {lower(obj.params.plotVarName)}, obj.params.maxVal, ...
                        obj.params.minVal, 'variableNames', {'expID', 'plotVarName', 'maxVal', ...
                        'minVal'});
            if isempty(obj.savedPlotParams)
                obj.savedPlotParams = newRow; 
            elseif isempty(obj.savedPlotParams(strcmp(obj.savedPlotParams.expID, newRow.expID) & ...
                        strcmpi(obj.savedPlotParams.plotVarName, newRow.plotVarName), :))
                obj.savedPlotParams = sortrows([obj.savedPlotParams; newRow]);
            else
                currRowIdx = strcmp(obj.savedPlotParams.expID, newRow.expID) & strcmpi( ...
                        obj.savedPlotParams.plotVarName, newRow.plotVarName);
                obj.savedPlotParams.maxVal(currRowIdx) = newRow.maxVal;
                obj.savedPlotParams.minVal(currRowIdx) = newRow.minVal;
            end
        end
        
        % PROCESS DATA AND CREATE THE SUMMARY PLOT
        function ax = plot_summary(obj, ax)
            
            % Determine whether this is a FicTrac or fluorescence summary
            if any(strcmpi(obj.params.plotVarName, {'moveSpeed', 'fwSpeed', 'yawSpeed', ...
                        'sideSpeed'}))
                dataType = 'ft';
            elseif any(strcmpi(obj.params.plotVarName, {'flailing', 'grooming', 'locomotion', ...
                        'isolatedmovement', 'ballstop', 'odor'}))
                dataType = 'annot';
            else
                dataType = 'fl';
            end
            
            params = obj.params; %#ok<*PROPLC>
            
            % Use saved colormap ranges if present
            if ~params.ignoreSavedParams && ~isempty(obj.savedPlotParams)
                savedParams = obj.savedPlotParams(strcmp(obj.savedPlotParams.expID, params.expID) & ...
                        strcmpi(obj.savedPlotParams.plotVarName, params.plotVarName), :);
                if ~isempty(savedParams)
                    params.maxVal = savedParams.maxVal;
                    params.minVal = savedParams.minVal;
                end
            end
            
            % Subset source DataTable appropriately
            filterDefs = struct();
            filterDefs.expID = params.expID;
            if strcmp(dataType, 'fl')
                filterDefs.roiName = params.plotVarName;
            else
                filterDefs.roiName = 'Background';
            end
            dt = obj.sourceDataTable.initialize_filters(filterDefs);
            
            % If there are bilateral ROIs, just use the first one
            uniqueROIs = unique(dt.subset.roiName);
            if numel(uniqueROIs) > 1
                filterDefs.roiName = uniqueROIs{1};
                dt = dt.initialize_filters(filterDefs);
                disp('Bilateral ROIs detected...using first ROI only')
            end
            dtSub = dt.subset;
            disp(unique(dtSub.roiName))
            
            % Specify sample rate
            if strcmp(dataType, 'fl')
                sampRate = dtSub.volumeRate(1);
            else
                sampRate = obj.FRAME_RATE;
            end
            
            % Extract and process data
            if strcmp(dtSub.groupName{1}, 'gaplessAcq')
                acqType = 'gapless';
            else
                acqType = 'singleTrial';
            end
            
            % Generate plotting data
            trialNumList = unique(dtSub.trialNum);
            blockedPlotData = {[]}; stimOnsets = {[]}; stimOffsets = {[]};
            for iTrial = 1:numel(trialNumList)
                currTrialData = dtSub(dtSub.trialNum == trialNumList(iTrial), :);
                
                % Extract the appropriate plot data
                if strcmp(dataType, 'fl')
                    currPlotData = currTrialData.rawFl{:};
                    if ~isempty(currTrialData.pmtShutoffVols{:}) & ...
                                ~isnan(currTrialData.pmtShutoffVols{:}) 
                            shutoffVols = currTrialData.pmtShutoffVols{:};
                            shutoffVols = shutoffVols(shutoffVols <= numel(currPlotData));
                            currPlotData(shutoffVols) = nan;
                    end
                    currPlotData = smoothdata(currPlotData, 1, 'gaussian', params.smWinVols, ...
                        'omitnan');
                else
                    if strcmp(dataType, 'annot')
                        currPlotData = zeros(sampRate * currTrialData.trialDuration, 1);
                    else
                        currPlotData = currTrialData.(params.plotVarName){:};
                        if strcmpi(params.plotVarName, 'yawSpeed')
                            currPlotData = rad2deg(currPlotData) * obj.FRAME_RATE;
                        else
                            currPlotData = currPlotData * 4.5 * obj.FRAME_RATE;
                        end
                        currPlotData = repeat_smooth(currPlotData, params.ftSmReps, 'dim', 1, ...
                                'smWin', params.smWinFrames);
                    end
                end
                
                % Limit data to max and min vals
                currPlotData(currPlotData > params.maxVal) = params.maxVal;
                currPlotData(currPlotData < params.minVal) = params.minVal;
                
                % Find stim onsets and offsets for current trial
                sampTimes = calc_volTimes(numel(currPlotData), sampRate, ...
                    currTrialData.trialDuration, currTrialData.originalTrialCount);
                eventNames = fieldnames(obj.eventData);
                allStimVec = zeros(size(sampTimes));
                for iName = 1:numel(eventNames)
                    currName = eventNames{iName};
                    if ismember(currName, params.stimEventNames)
                        eventArr = obj.eventData.(currName).create_logical_array( ...
                            numel(sampTimes), sampTimes, currTrialData(:, {'expID', 'trialNum'}));
                        allStimVec = allStimVec | eventArr';
                    end
                end
                currStimOnsets = zeros(size(sampTimes));
                currStimOffsets = zeros(size(sampTimes));
                currStimOnsets(strfind(allStimVec, [0 1])) = 1;
                currStimOffsets(strfind(allStimVec, [1 0])) = 1;
                if allStimVec(end) == 1
                    currStimOffsets(end) = 1;
                end
                if allStimVec(1) == 1
                    currStimOnsets(1) = 1;
                end
                
                % Generate logical array for flailing if necessary
                if strcmp(dataType, 'annot')
                    annotEventName = obj.params.plotVarName;
                    annotEventVec = obj.eventData.(annotEventName).create_logical_array( ...
                        numel(sampTimes), sampTimes, currTrialData(:, {'expID', 'trialNum'}));
                    currPlotData = annotEventVec;
                end
                
                % Save plot and stim timing data as one or more blocks, reshaping as necessary
                currPlotData = reshape(currPlotData, [], currTrialData.originalTrialCount);
                currStimOnsets = reshape(currStimOnsets, [], currTrialData.originalTrialCount);
                currStimOffsets = reshape(currStimOffsets, [], currTrialData.originalTrialCount);
                if strcmp(acqType, 'gapless')
                    blockedPlotData{iTrial} = currPlotData;
                    stimOnsets{iTrial} = currStimOnsets;
                    stimOffsets{iTrial} = currStimOffsets;
                else
                    blockedPlotData{1} = [blockedPlotData{1}, currPlotData];
                    stimOnsets{1} = [stimOnsets{1}, currStimOnsets];
                    stimOffsets{1} = [stimOffsets{1}, currStimOffsets];
                end
                                
            end%iTrial
            
            % Plot data
            allPlotData = cell2mat(blockedPlotData)';
            trialGroups = repelem(1:numel(blockedPlotData), cellfun(@(x) size(x, 2), ...
                    blockedPlotData));
            plot_2D_summary(allPlotData, sampRate, 'plotaxes', ax, 'trialGroups', ...
                    trialGroups, 'spacerArrRows', params.spacerRows);
            colorbar();
            
            
            % Add overlay of stim timing lines
            hold on;
            yCount = 0.5;
            yL = ylim();
            for iBlock = 1:numel(stimOnsets)
                currOnsets = stimOnsets{iBlock};
                currOffsets = stimOffsets{iBlock};
                for iTrial = 1:size(currOnsets, 2)
                    currOnsetSamps = find(currOnsets(:, iTrial));
                    currOffsetSamps = find(currOffsets(:, iTrial));
                    yy = [yCount, yCount + 1];
                    for iStim = 1:numel(currOnsetSamps)
                        xx = [1 1] * currOnsetSamps(iStim);
                        plot(xx, yy, 'linewidth', 2, 'color', 'b');
                        xx = [1 1] * currOffsetSamps(iStim);
                        plot(xx, yy, 'linewidth', 2, 'color', 'r');
                    end
                    yCount = yCount + 1;
                end
                yCount = yCount + params.spacerRows;
            end
            ylim(yL);
            titleStr = params.plotVarName;
            if strcmp(dataType, 'fl')
                titleStr = ['Raw F — ', titleStr];
            end
            title(titleStr);
        end
        
    end%Methods
    
end%Class

% ==================================================================================================
% Local functions
% ==================================================================================================
function paramStruct = initialize_default_params(obj)
    paramStruct = struct();
    paramStruct.expID = obj.expList.expID{1};
    paramStruct.plotVarName = 'TypeD';
    paramStruct.stimEventNames = {'odor', 'soundstim', 'optostim'};
    paramStruct.maxVal = 1500;
    paramStruct.minVal = 0;
    paramStruct.smWinVols = 3;
    paramStruct.smWinFrames = 3;
    paramStruct.ftSmReps = 10;
    paramStruct.spacerRows = 2;
    paramStruct.ignoreSavedParams = 0;
end
