classdef behavAnnotEvent < event
   
    properties (SetAccess = immutable)
        metadataFieldNames
    end
    
    methods
        
        % Constructor
        function obj = behavAnnotEvent(expID, type)
            if ~ismember(type, {'locomotion', 'isolatedMovement', 'grooming', 'quiescence'})
               error('invalid annotation type'); 
            end
            obj = obj@event(expID, type);
        end
        
        function obj = append_annot_data(obj, trialNums, annotData, frameTimes)
           
            % Make sure input variable sizes match up
            if numel(trialNums > 1)
                if isnumeric(annotData) || isnumeric(frameTimes)
                    error('when appending multiple trials annotData must be in a cell array')
                elseif numel(trialNums) ~= numel(annotData) || numel(trialNums) ~= numel(frameTimes)
                    error('size mismatch between trialNums and annotData or frameTimes')
                end
            elseif numel(trialNums) == 1 && ~iscell(annotData)
                annotData = {annotData}; % Put flow data in a cell if necessary
            elseif numel(trialNums) == 1 && ~iscell(frameTimes)
                frameTimes = {frameTimes};
            end
            
            for iTrial = 1:numel(trialNums)
                
                currAnnotData = annotData{iTrial};                
                if sum(currAnnotData) > 0
                    % Identify event onset frames
                    onsetFrames = find(diff(currAnnotData) == 1) + 1;
                    if currAnnotData(1) == 1
                        onsetFrames = [1, onsetFrames];
                    end
                    
                    % Identify move offset frames
                    offsetFrames = find(diff(currAnnotData) == -1) + 1;
                    if currAnnotData(end) == 1
                        offsetFrames = [offsetFrames, numel(currAnnotData)];
                    end
                    
                    % Convert to times
                    onsetTimes = frameTimes{iTrial}(onsetFrames);
                    offsetTimes = frameTimes{iTrial}(offsetFrames);
                    eventTimes = {[onsetTimes', offsetTimes']};
                    
                    mdFieldNames = {'moveThresh', 'smReps', 'smWin'};
                    mdFieldVals = {threshold, 20, 6};
                    
                    % Append data
                    obj = obj.append_data(trialNums(iTrial), eventTimes, mdFieldNames, mdFieldVals);
                end
            end%iTrial
            
            
            
        end
        
        
    end
    
end