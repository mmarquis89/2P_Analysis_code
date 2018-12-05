function calc_event_dff_avg(parentDir, filterStr, sessionDataFile)
%===================================================================================================
% Calculate volume- and trial-averaged dF/F for a set of event data
%
% Required inputs:
%       A filter string using * as a wildcard that identifies the files in parentDir to be processed.
%       The name of an imaging data file in the parent directory.
%       That file should contain the following variables:
%           alignEventSummary
%           filterEventSummary
%           primaryEventNames
%           eventLists
%           nEventTypes
%           condNames
%           onsetFilterVecs
%           offsetFilterVecs
%           
%       Output file names are based on the input file names.
%
% Outputs:
%       A .mat file for each input file containing the following variables:
%           combinedDffAvg
%           allCondSummaries
%           allCondNames
%           combFilterVecs
%===================================================================================================
try
addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

disp(parentDir)
disp(sessionDataFile)

% Identify event data files
eventDataFiles = dir(fullfile(parentDir, filterStr));

% Load session data file
[analysisMetadata, wholeSession] = load_imaging_data(parentDir, sessionDataFile, 'LoadSessionData', 1);
sessionSize = size(wholeSession);

for iFile = 1:numel(eventDataFiles)
    
    % Load event data variables
    load(fullfile(parentDir, eventDataFiles(iFile).name));
    disp(eventDataFiles(iFile).name);
    whos
%           alignEventSummary
%           filterEventSummary
%           primaryEventNames
%           eventLists
%           nEventTypes
%           condNames
%           onsetFilterVecs
%           offsetFilterVecs
%           analysisWindows
    

    % ----------------------------------------------------------------------------------------------
    % Calculate dF/F
    % ----------------------------------------------------------------------------------------------
    
    onsetDffAvg = []; offsetDffAvg = []; onsetDff = []; offsetDff = []; combinedDff = []; combinedDffAvg =[]; respRawF = []; baselineRawF = [];
    onsetCondSummaries = []; offsetCondSummaries = []; allCondSummaries = [];
    for iType = 1:nEventTypes
        
        primaryFiltName = primaryEventNames{iType};
        analysisWindow = analysisWindows(iType, :);
        nConds = numel(condNames{iType});
        eventList = eventLists{iType};
        
        baselineDur = analysisWindow(1);
        respDur = analysisWindow(2);
        onsetDffAvg{iType} = zeros([sessionSize(1:3), nConds]);
        offsetDffAvg{iType} = onsetDffAvg{iType};
        
        for iCond = 1:nConds
            
            disp(['Calculating dF/F for ', primaryFiltName, ' cond #', num2str(iCond), ' of ', num2str(nConds), '...'])
            
            % Calculate dF/F for event onsets
            if sum(onsetFilterVecs{iType}(:,iCond)) > 0
                
                [baselineData, respData] = extract_event_volumes(eventList, onsetFilterVecs{iType}(:,iCond), baselineDur, respDur, analysisMetadata, ...
                    wholeSession, 'offsetAlign', 0); % --> [y, x, plane, volume, event]
                
                currDffAvg = calc_dFF(respData, baselineData, [4 5]);                    % --> [y, x, plane]
                
                baselineRawF{iType}{iCond} = baselineData;                               % --> {eventType}{Condition}[y, x, plane, volume, event]
                respRawF{iType}{iCond} = respData;                                       % --> {eventType}{Condition}[y, x, plane, volume, event]
                
                onsetDffAvg{iType}(:,:,:, iCond) = currDffAvg;                           % --> {eventType}[y, x, plane, condition]
                
            end
            
            % Calculate dF/F for event offsets
            if sum(offsetFilterVecs{iType}(:,iCond)) > 0
                
                [baselineData, respData] = extract_event_volumes(eventList, offsetFilterVecs{iType}(:,iCond), baselineDur, respDur, analysisMetadata, ...
                    wholeSession, 'offsetAlign', 1); % --> [y, x, plane, volume, event]
                
                currDffAvg = calc_dFF(respData, baselineData, [4 5]);                       % --> [y, x, plane]
                
                offsetDffAvg{iType}(:,:,:, iCond) = currDffAvg;                             % --> {eventType}[y, x, plane, condition]
            end
        end% iCond
        clear respData baselineData
        
        combinedDffAvg{iType} = cat(4, onsetDffAvg{iType}, offsetDffAvg{iType});            % --> {eventType}[y, x, plane, condition]
        clear  onsetDffAvg offsetDffAvg
        
        % Create summary table for onset conditions
        condCountCol = (sum(onsetFilterVecs{iType})');
        alignCol = repmat({'onset'}, nConds, 1);
        baselineCol = repmat(sprintf('%g', analysisWindow(1)), nConds, 1);
        respCol = repmat(sprintf('%g', analysisWindow(2)), nConds, 1);
        rowNames = cellfun(@num2str, num2cell(1:nConds), 'uniformOutput', 0);
        varNames = {'Count', 'CondName', 'Align', 'Base', 'Resp'};
        onsetCondSummaries{iType} = table(condCountCol, condNames{iType}, alignCol, baselineCol, respCol, 'RowNames', rowNames, 'VariableNames', varNames);
        
        % Create summary table for offset conditions
        condCountCol = (sum(offsetFilterVecs{iType})');
        alignCol = repmat({'offset'}, nConds, 1);
        baselineCol = repmat(sprintf('%g', analysisWindow(1)), nConds, 1);
        respCol = repmat(sprintf('%g', analysisWindow(2)), nConds, 1);
        rowNames = cellfun(@num2str, num2cell((1:nConds) + nConds), 'uniformOutput', 0);
        varNames = {'Count', 'CondName', 'Align', 'Base', 'Resp'};
        offsetCondSummaries{iType} = table(condCountCol, condNames{iType}, alignCol, baselineCol, respCol, 'RowNames', rowNames, 'VariableNames', varNames);
        
        % Concatenate onset and offset summary tables
        allCondSummaries{iType} = vertcat(onsetCondSummaries{iType}, offsetCondSummaries{iType});
        allCondNames{iType} = [condNames{iType}; condNames{iType}];
        
        % Remove conditions with one or fewer occurences
        nullConds = allCondSummaries{iType}.Count <= 1; % using one instead of zero because it causes a bug later on
        combFilterVecs{iType} = [onsetFilterVecs{iType}, offsetFilterVecs{iType}];
        combFilterVecs{iType}(:,nullConds) = [];
        allCondSummaries{iType}(nullConds, :) = [];
        allCondSummaries{iType}.Properties.RowNames = cellfun(@num2str, num2cell(1:size(allCondSummaries{iType}, 1)), 'UniformOutput', 0);
        allCondNames{iType}(nullConds) = [];
        combinedDffAvg{iType}(:,:,:, nullConds) = [];   % --> {eventType}[y, x, plane, condition]
        disp(allCondSummaries{iType})
        
    end% iType
    disp('dF/F calculation complete')
    
    % Save combinedDffAvg and summaries
    save(fullfile(parentDir, ['CombDffAvg_', eventDataFiles(iFile).name]), 'combinedDffAvg', 'allCondSummaries', ...
        'allCondNames', 'combFilterVecs', '-v7.3')
    
end% iFile
catch ME
    write_to_log(ME.message, mfilename);
end%try
end% function