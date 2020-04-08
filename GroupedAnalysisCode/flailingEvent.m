classdef flailingEvent < event

    
    properties (SetAccess = immutable)
        metadataFieldNames
    end
    
    methods
        
        % Constructor
        function obj = flailingEvent(expID)
            type = 'flailing';
            obj = obj@event(expID, type);
        end
        
        % Modified append_data function for flailing
        function obj = append_flow_data(obj, trialNums, flowData, frameTimes, threshold)
            
            % Make sure input variable sizes match up
            if numel(trialNums > 1)
                if isnumeric(flowData) || isnumeric(frameTimes)
                    error('when appending multiple trials flow data must be in a cell array')
                elseif numel(trialNums) ~= numel(flowData) || numel(trialNums) ~= numel(frameTimes)
                    error('size mismatch between trialNums and flowData or frameTimes')
                end
            elseif numel(trialNums) == 1 && ~iscell(flowData)
                flowData = {flowData}; % Put flow data in a cell if necessary 
            elseif numel(trialNums) == 1 && ~iscell(frameTimes)
                frameTimes = {frameTimes};
            end
            
            for iTrial = 1:numel(trialNums)
                
                % Process and threshold flow data
                smFlow = repeat_smooth(flowData{iTrial}, 20, 'dim', 2, 'smwin', 6);
                smFlow = smFlow - min(smFlow);
                moveFrames = smFlow > threshold;
                
                if sum(moveFrames) > 0
                    % Identify move onset frames
                    moveOnsetFrames = find(diff(moveFrames) == 1) + 1;
                    if moveFrames(1) == 1
                        moveOnsetFrames = [1, moveOnsetFrames];
                    end
                    
                    % Identify move offset frames
                    moveOffsetFrames = find(diff(moveFrames) == -1) + 1;
                    if moveFrames(end) == 1
                        moveOffsetFrames = [moveOffsetFrames, numel(moveFrames)];
                    end
                    
                    % Convert to times
                    moveOnsetTimes = frameTimes{iTrial}(moveOnsetFrames);
                    moveOffsetTimes = frameTimes{iTrial}(moveOffsetFrames);
                    eventTimes = {[moveOnsetTimes', moveOffsetTimes']};
                    
                    mdFieldNames = {'moveThresh', 'smReps', 'smWin'};
                    mdFieldVals = {threshold, 20, 6};
                    
                    % Append data
                    obj = obj.append_data(trialNums(iTrial), eventTimes, mdFieldNames, mdFieldVals);
                end
            end%iTrial

        end%function
        
        
        function obj = append_annotation_data(obj)
            
            
        end
        
        
    end%methods    
    
end