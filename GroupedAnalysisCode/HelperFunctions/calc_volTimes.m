function volTimes = calc_volTimes(nVolumes, volumeRate, trialDuration, originalTrialCount)

% Special volTimes calculation for compatibility with Gapless acq experiments. Works with other 
% experiments too since the "originalTrialCount" will equal one, and therefore the calculation will 
% simplify to: (1:nVolumes) / volumeRate.

volsPerTrial = nVolumes / originalTrialCount;
trialVolTimes = (1:volsPerTrial) / volumeRate;
volTimes = [];
for iTrial = 1:originalTrialCount
    volTimes = [volTimes, (trialVolTimes + ((trialDuration / ...
        originalTrialCount) * (iTrial - 1)))];
end


end