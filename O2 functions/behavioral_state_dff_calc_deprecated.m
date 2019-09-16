function behavioral_state_dff_calc(parentDir, sessionDataFile)

% CALCULATE AND PLOT OVERALL MEAN dF/F ACROSS BEHAVIORAL STATES

try
addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

% Load data
[analysisMetadata, wholeSession] = load_imaging_data(parentDir, sessionDataFile, 'LoadSessionData', 1);

behaviorNames = {'Locomotion','IsoMove', 'AnyMove'};
behaviorNums = {3, 1, [1 3]};
baselineLabel = 0;

for iType = 1:numel(behaviorNames)
    
    actionLabel = behaviorNums{iType};
    
    % Identify behavioral state during each volume
    actionVols = zeros(analysisMetadata.nTrials, analysisMetadata.nVolumes); stoppedVols = actionVols;
    for iTrial = 1:analysisMetadata.nTrials
        if analysisMetadata.goodTrials(iTrial)
            
            % Pull out action numbers for each volume
            currActions = analysisMetadata.trialAnnotations(iTrial, :);
            volActions = currActions(analysisMetadata.volFrames);
            
            % Identify volume actions
            actionVols(iTrial, :) =  ismember(volActions, actionLabel);   %--> [trial, vol]
            stoppedVols(iTrial, :) = ismember(volActions, baselineLabel); %--> [trial, vol]
        else
            % So data from invalid trials won't ever be matched to an action state
            actionVols(iTrial, :) = 0;
            stoppedVols(iTrial, :) = 0;
        end
    end
    
    % Calculate average values for each plane across behavioral states
    meanActionVols = [];
    meanStoppedVols = [];
    for iTrial = 1:analysisMetadata.nTrials
        disp(['Trial ', num2str(iTrial)])
        currImgData = wholeSession(:,:,:,:,iTrial);
        currActionVols = logical(actionVols(iTrial,:));
        currStoppedVols = logical(stoppedVols(iTrial,:));
        
        % Pull out running volumes, if any exist, from the current trial
        if sum(currActionVols) > 0
            meanActionVols(:,:,:,end+1) = mean(currImgData(:,:,:,currActionVols),4);    %--> [y, x, plane, trial]
        end
        
        % Pull out stopping volumes, if any exist, from the current trial
        if sum(currStoppedVols) > 0
            meanStoppedVols(:,:,:,end+1) = mean(currImgData(:,:,:,currStoppedVols),4);  %--> [y, x, plane, trial]
        end 
    end
    
    actionMean = mean(meanActionVols, 4);          %--> [y, x, plane]
    stoppedMean = mean(meanStoppedVols, 4);        %--> [y, x, plane]
    
    % Get dF/F values for action relative to quiescence
    actionDff = (actionMean - stoppedMean) ./ stoppedMean; % --> [y, x, plane]
    
    % Save dFF data
    save(fullfile(parentDir, ['actionDff_', behaviorNames{iType}]), 'actionDff', '-v7.3')
end
catch ME
    write_to_log(getReport(ME), mfilename);
end%try
end