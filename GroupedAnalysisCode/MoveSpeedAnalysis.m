classdef MoveSpeedAnalysis
% ==================================================================================================   
%    
% Properties:
%
%       sourceData
%       eventObjects
%       filterDefs
%       analysisParams
%       
% Methods:
%
% 
%
% ==================================================================================================    
    
properties
    sourceData % DataTable with all data except events
    eventObjects
    filterDefs
    analysisParams
end

properties (Dependent)
    eventNames
end

methods
    
    % Constructor
    function obj = MoveSpeedAnalysis(expList)
        
        % Load metadata, FicTrac data, ROI data for all experiments
        parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';
        allExpData = [];
        for iExp = 1:size(expList, 1)
            
            currExpID = expList.expID{iExp};
            disp(currExpID);
            
            % Generate full paths to data files
            roiFileName = [currExpID, '_roiData.mat'];
            roiFile = fullfile(parentDir, roiFileName);
            ftFileName = [currExpID, '_ficTracData.mat'];
            ftFile = fullfile(parentDir, ftFileName);
            expMdFile = fullfile(parentDir, [currExpID, '_expMetadata.csv']);
            trialMdFile = fullfile(parentDir, [currExpID, '_trialMetadata.mat']);
            
            if exist(roiFile, 'file') && exist(ftFile, 'file') && exist(expMdFile, 'file')
                load(roiFile, 'roiData');
                load(ftFile, 'ftData');
                expMd = readtable(expMdFile, 'delimiter', ',');
                load(trialMdFile, 'trialMetadata');
                trialMd = trialMetadata;
                if exist(odorEventFile, 'file')
                    odorEvents = readtable(odorEventFile, 'delimiter', ',');
                    currData = inner_join(expMd, trialMd, roiData, ftData);
                    allExpData = [allExpData; outerjoin(currData, odorEvents, 'type', 'left', ...
                        'mergekeys', 1)];
                else
                    allExpData = [allExpData; inner_join(expMd, trialMd, roiData, ftData)];
                    
                end
            else
                disp(['Skipping ', currExpID, ' due to one or more missing files']);
            end
        end%iExp
        
        % Create source DataTable
        obj.sourceData = DataTable(allExpData);
        obj.filterDefs = struct();
        obj.filterDefs.roiName = [];
        obj.filterDefs.expID = [];
        obj.sourceData = obj.sourceData.initialize_filters(obj.filterDefs);
        
        % Load all event data objects
        obj.eventObjects = load_event_data(expList, parentDir);
        
        % Set initial analysis parameters
        obj.analysisParams = initialize_analysis_params();
    end
    
    % Dependent property methods
    function eventNames = get.eventNames(obj)
        eventNames = fieldNames(obj.eventObjects)';
    end
    
    % Set target ROI name and expID
    function obj = set_exp_filters()
    
    end
    
    
    
    
end%methods

end%class

% ==================================================================================================
% Local functions
% ==================================================================================================
function paramStruct = initialize_analysis_params(obj)
    paramStruct = struct();
    paramStruct.max_moveSpeed = 100;
    paramStruct.smWinVols = 1;
    paramStruct.smWinFrames = 1;
    paramStruct.nSmoothRepsFrames = 15;
    paramStruct.lagVols = 2;
    paramStruct.nHistBins = 50;
    paramStruct.plotting_minSpeed = 0;
    
    paramStruct.excludeEvents = struct();
    for iType = 1:numel(obj.eventNames)
        paramStruct.(obj.eventNames{iType}) = 0;
    end
end













