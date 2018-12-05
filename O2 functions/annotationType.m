classdef annotationType
    properties
        % From construction
        frameAnnotArr
        volAnnotArr
        linAnnotArr
        linAnnotArrStr
        name
        
        % From methods
        onsetVolsLin
        offsetVolsLin
        eventVolsLin
        onsetVols
        offsetVols
        eventVols
        eventList
        nEvents
        
    end
    methods
        % Constructor function
        function obj = annotationType(infoStruct, annotArr, name)
            if nargin > 0
                obj.name = name;
                obj.frameAnnotArr = annotArr;
                    obj.frameAnnotArr(:, [1 size(obj.frameAnnotArr, 2)]) = 0;
                for iTrial = 1:infoStruct.nTrials
                    obj.volAnnotArr(iTrial, :) = obj.frameAnnotArr(iTrial, infoStruct.volFrames');
                end
                    obj.volAnnotArr(:, [1 size(obj.volAnnotArr, 2)]) = 0;
                obj.linAnnotArr = annot2lin(obj.volAnnotArr);
                obj.linAnnotArrStr = num2str(obj.linAnnotArr);
                    obj.linAnnotArrStr = obj.linAnnotArrStr(~isspace(obj.linAnnotArrStr));
            end
        end
        
        % Get onset, offset, and event volumes, then extract list
        function obj = get_event_vols(obj, onsetStr, offsetStr)
            
            % Determine onset and offset frames (linearized) for all events and convert to logical vectors
            obj.onsetVolsLin = zeros(size(obj.linAnnotArr));
            obj.onsetVolsLin(regexp(obj.linAnnotArrStr, onsetStr) + 1) = 1;
            obj.offsetVolsLin = zeros(size(obj.linAnnotArr));
            obj.offsetVolsLin(regexp(obj.linAnnotArrStr, offsetStr) + 1) = 1;
            obj.eventVolsLin = zeros(size(obj.linAnnotArr));
            obj.eventVolsLin(cell2mat(arrayfun(@(x, y) {x:y}, regexp(obj.linAnnotArrStr, onsetStr), ...
                                      regexp(obj.linAnnotArrStr, offsetStr)))) = 1; 
                                      
            % Convert logical vectors back into 2D arrays                          
            arrSize = size(obj.volAnnotArr);
            obj.onsetVols = logical(lin2annot(obj.onsetVolsLin, arrSize));
            obj.offsetVols = logical(lin2annot(obj.offsetVolsLin, arrSize));
            obj.eventVols = logical(lin2annot(obj.eventVolsLin, arrSize));
            
            % Create chronological lists of events
            obj.eventList = create_event_list(obj.onsetVols, obj.offsetVols);
            obj.nEvents = size(obj.eventList, 1);
        end       
        
    end    
end