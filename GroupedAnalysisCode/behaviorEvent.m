classdef behaviorEvent < event
    
    properties
        metadataFieldNames
    end
    
    methods
        % Constructor
        function obj = behaviorEvent(expID, type)
            if ~ismember(type, {'locomotion', 'isolatedMovement', 'grooming', 'quiescence', ...
                        'flailing'})
                error('invalid event type');
            end
            obj = obj@event(expID, type);
            obj.metadataFieldNames = [];
        end
        
        function obj = append_annotation_data(obj, trialNums, annotData, frameTimes, varargin)
            
            % Parse optional arguments
            p = inputParser;
            addParameter(p, 'MetadataFieldNames', []);
            addParameter(p, 'MetadataFieldValues', []);
            parse(p, varargin{:});
            mdFieldNames = p.Results.MetadataFieldNames;
            mdFieldVals = p.Results.MetadataFieldValues;
            
            % Make sure input variable sizes and types match up
            [annotData, frameTimes] = multi_trial_input_check(trialNums, annotData, frameTimes);
            
            % Check for match with metaDataFieldNames
            if ~isempty(obj.metadataFieldNames)
                if ~strcmp(obj.metadataFieldNames, mdFieldNames)
                    error('metadata fields are incompatible with previously existing data')
                end
            else
                obj.metadataFieldNames = mdFieldNames;
            end           
            
            % Identify event times and append eventData
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
                    
                    % Append data
                    obj = obj.append_data(trialNums(iTrial), eventTimes, mdFieldNames, mdFieldVals);
                end
                
            end%iTrial
        end%function
            
        
        function obj = append_flow_data(obj, trialNums, flowData, frameTimes, threshold)           
            
            % Make sure input variable sizes and types match up
            [flowData, frameTimes] = multi_trial_input_check(trialNums, flowData, frameTimes);
            
            % Threshold flow data
            for iTrial = 1:numel(trialNums)

                smFlow = repeat_smooth(flowData{iTrial}, 20, 'dim', 2, 'smwin', 6);
                smFlow = smFlow - min(smFlow);
                moveFrames = smFlow > threshold;
                
                mdFieldNames = {'moveThresh', 'smReps', 'smWin'};
                mdFieldVals = {threshold, 20, 6};
                
                if sum(moveFrames) > 0
                    obj = obj.append_annotation_data(trialNums(iTrial), moveFrames, ...
                            frameTimes{iTrial}, 'MetadataFieldNames', mdFieldNames, ...
                            'MetadataFieldValues', mdFieldVals);
                end 
            end
        end%function
        
    end%methods
end %class


% ==================================================================================================
% Local functions
% ==================================================================================================


function [behavData, frameTimes] = multi_trial_input_check(trialNums, behavData, frameTimes) 

    % Make sure input variable sizes and types match up
    if numel(trialNums) > 1
        if isnumeric(behavData) || isnumeric(frameTimes)
            error(['when appending multiple trials behavior data and frameTimes must be in cell', ...
                    ' arrays'])
        elseif numel(trialNums) ~= numel(behavData) || numel(trialNums) ~= numel(frameTimes)
            error('size mismatch between trialNums and the behavior data or frameTimes')
        end
    else
        if ~iscell(behavData)
            behavData = {behavData}; % Put flow data in a cell if necessary
        end
        if ~iscell(frameTimes)
            frameTimes = {frameTimes};
        end
    end
end 













