classdef RegressionModelAnalysis_PPM1201 < RegressionModelAnalysis
    
    methods
        
        % Constructor
        function obj = RegressionModelAnalysis_PPM1201(expInfoTbl, sourceDataParams)
            
        % expInfoTable = table with colums: [expID][skipTrials][skipVols]
        
        % sourceDataParams = struct with fields: roiName 
        %                                        maxSpeed 
        %                                        smWinVols
        %                                        smWinFrames
        %                                        smReps 
        %                                        ftLagVols 
        %                                        speedType 
        %                                        odorRespFilterDur
        %                                        parentDir
        %                                        dataDir
        %                                        eventDataParentDir
        %                                        alignEventDateStr
        %                                        convertFtUnits
        %                                        loadOneNoteData
        %                                        alignObjFilterDefs            
            p = sourceDataParams;
            
            % Set default value for loading oneNoteData
            if ~isfield(p.loadOneNoteData) || isempty(p.loadOneNoteData)
                p.loadOneNoteData = 0;
            end
            
            % Set alignObj filterDefs manually
            if ~isfield(p.alignObjFilterDefs) || isempty(p.alignObjFilterDefs)
                p.alignObjFilterDefs = ...
                        alignObj.create_filterDefs('loadOneNoteData', p.loadOneNoteData);
                 p.alignObjFilterDefs.expID = '20180309-2';
                 p.alignObjFilterDefs.roiName = 'TypeF-R';
                 p.alignObjFilterDefs.odor = -1;
                 p.alignObjFilterDefs.locomotion = 0;
            end
            
            % Make sure FicTrac data units are converted
            if  ~isfield(p.convertFtUnits) || isempty(p.convertFtUnits)
                p.convertFtUnits = 1;
            end
            
            % Initialize object
            obj = obj@RegressionModelAnalysis(expInfoTbl, p);            
            
        end        
    end%Methods
    
end%Class