function [annotationTypes, annotationTypeSummary] = process_annotation_types(analysisMetadata)

try
    
% Copy variables for convenience
nVolumes = analysisMetadata.nVolumes;
if ~isempty(analysisMetadata.nFrames)
    nFrames = analysisMetadata.nFrames;
else
    nFrames = nVolumes;
end
nTrials = analysisMetadata.nTrials;
write_to_log(num2str(nTrials), mfilename);
write_to_log(num2str(size(analysisMetadata.goodTrials)));
goodTrials = analysisMetadata.goodTrials(1:nTrials);

%   annotArr: (row = trial, col = frame, Z-dim = event type)
%   Event types are: [open loop stim, behavior]

annotationTypes = [];

% ----------------------------------------------------------------------------------------------
% Closed loop odor stim
% ----------------------------------------------------------------------------------------------
if contains(analysisMetadata.stimTypes, 'ClosedLoop')
    write_to_log(num2str(size(analysisMetadata.daqFtData)), mfilename);
    stimAnnotArrCL = round(squeeze(analysisMetadata.daqFtData(5,:,:))'); % --> [trial, sample] Copy of the stim command signal
    write_to_log(num2str(size(stimAnnotArrCL)), mfilename);
    odorAnnotationsCL = annotationType(analysisMetadata, stimAnnotArrCL, 'ClosedLoopStim');
    odorAnnotationsCL = get_event_vols(odorAnnotationsCL, '05', '50');
    annotationTypes{end + 1} = odorAnnotationsCL;
end

% ----------------------------------------------------------------------------------------------
% Open loop stimuli
% ----------------------------------------------------------------------------------------------

% Add stim frames
stimTypes = analysisMetadata.stimTypes(~contains(analysisMetadata.stimTypes, 'ClosedLoop'));
nStims = numel(stimTypes);
annotArr = zeros(nTrials, nFrames, nStims);
for iStim = 1:nStims
    for iTrial = 1:nTrials
        onsetFrame = round(analysisMetadata.stimOnsetTimes(iTrial)) * analysisMetadata.FRAME_RATE;
        if contains(analysisMetadata.stimTypes{iStim}, 'Pair')
            % Extend event window to cover second stimulus
            offsetFrame = (round(analysisMetadata.stimOnsetTimes(iTrial)) + round(3 * analysisMetadata.stimDurs(iTrial))) * analysisMetadata.FRAME_RATE;
        else
            offsetFrame = (round(analysisMetadata.stimOnsetTimes(iTrial)) + round(analysisMetadata.stimDurs(iTrial))) * analysisMetadata.FRAME_RATE;
        end
        if analysisMetadata.stimSepTrials.(stimTypes{iStim})(iTrial)
            annotArr(iTrial, onsetFrame:offsetFrame, iStim) = 4; % --> [trial, frame]
        end
    end
end
annotArr(~goodTrials, :, :) = 0;

% Create annotation type for each stim
for iStim = 1:nStims
    currAnnotArr = annotArr(:, :, iStim);
    currAnnotations = annotationType(analysisMetadata, currAnnotArr, stimTypes{iStim});
    currAnnotations = get_event_vols(currAnnotations, '04', '40');
    annotationTypes{end + 1} = currAnnotations;
end

% ----------------------------------------------------------------------------------------------
% Behavior
% ----------------------------------------------------------------------------------------------

% Add behavior annotations
behaviorAnnotArr = analysisMetadata.trialAnnotations; % --> [trial, frame]
behaviorAnnotArr(~goodTrials, :) = 0;

%--------------------------------------------------------
% Quiescence = 0, Isolated movement = 1, Locomotion = 3
%--------------------------------------------------------

% Locomotion events
locomotionAnnotations = annotationType(analysisMetadata, behaviorAnnotArr, 'locomotion');
locomotionAnnotations = get_event_vols(locomotionAnnotations, '[01]3', '3[01]');
annotationTypes{end + 1} = locomotionAnnotations;

% Isolated movement events
isoMoveAnnotations = annotationType(analysisMetadata, behaviorAnnotArr, 'isoMove');
isoMoveAnnotations = get_event_vols(isoMoveAnnotations, '[03]1', '1[03]');
annotationTypes{end + 1} = isoMoveAnnotations;

% All behavior events
onsetRegExpStr = '(?<=[03])1|(?<=[01])3';
offsetRegExpStr = '1(?=[03])|3(?=[01])';
disp(['behaviorAnnotArr: ', num2str(size(behaviorAnnotArr))])
disp(onsetRegExpStr)
disp(offsetRegExpStr)
behaviorAnnotations = annotationType(analysisMetadata, behaviorAnnotArr, 'move');
behaviorAnnotations = get_event_vols(behaviorAnnotations, onsetRegExpStr, offsetRegExpStr);
annotationTypes{end + 1} = behaviorAnnotations;

% ==================================================================================================

% Make list of all annotation type names
for iType = 1:numel(annotationTypes)
    annotationTypeNames{iType} = annotationTypes{iType}.name;
end
annotationTypeSummary = table((1:numel(annotationTypeNames))', annotationTypeNames', 'VariableNames', {'Index', 'AnnotationType'});

catch ME
    write_to_log(getReport(ME), mfilename);
end%try
end