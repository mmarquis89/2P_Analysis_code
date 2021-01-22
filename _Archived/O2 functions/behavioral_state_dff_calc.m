function behavioral_state_dff_calc(parentDir, sessionDataFile)

% CALCULATE AND PLOT OVERALL MEAN dF/F ACROSS BEHAVIORAL STATES

try
addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

% Load data
[mData, ~] = load_imaging_data(parentDir, sessionDataFile);
nVolumes = mData.nVolumes;
nTrials = mData.nTrials;

behaviorNames = {'Locomotion','IsoMove', 'AnyMove'};
behaviorNums = {3, 1, [1 3]};
baselineLabel = 0;

actionDff = [];
for iPlane = 1:mData.nPlanes
    load(fullfile(parentDir, ['rigid_reg_sid_', num2str(mData.sid), ...
            '_chan_1_plane_', num2str(iPlane), '_sessionFile.mat'])); % --> 'sessionData'

    for iType = 1:numel(behaviorNames)

        actionLabel = behaviorNums{iType};

        % Identify behavioral state during each volume
        actionVols = zeros(nTrials, nVolumes); stoppedVols = actionVols;
        for iTrial = 1:nTrials
            if mData.goodTrials(iTrial)

                % Pull out action numbers for each volume
                currActions = mData.trialAnnotations(iTrial, :);
                volActions = currActions(mData.volFrames);

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
        actionMean = mean(sessionData(:,:, logical(actionVols(:))), 3, 'omitnan'); % --> [y, x]
        stoppedMean = mean(sessionData(:,:, logical(stoppedVols(:))), 3, 'omitnan'); % --> [y, x]

        % Get dF/F values for action relative to quiescence
        actionDff(:, :, iPlane, iType) = (actionMean - stoppedMean) ./ stoppedMean; % --> [y, x, plane]
    end%iType
end%iPlane

        % Save dFF data
        save(fullfile(parentDir, 'actionDff.mat'), 'actionDff', 'behaviorNames', '-v7.3')
        
catch ME
    write_to_log(getReport(ME), mfilename);
end%try
end