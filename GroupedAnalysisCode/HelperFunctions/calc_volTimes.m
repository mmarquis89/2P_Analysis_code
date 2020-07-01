function sampleTimes = calc_volTimes(nSamples, sampleRate, trialDuration, originalTrialCount)
% sampleTimes = calc_volTimes(nSamples, sampleRate, trialDuration, originalTrialCount)
% Special volTimes/frameTimes calculation for compatibility with Gapless acq experiments. Works with
% other experiments too since the "originalTrialCount" will equal one, and therefore the calculation 
% will simplify to: (1:nSamples) / sampleRate.

samplesPerTrial = nSamples / originalTrialCount;
trialSampleTimes = (1:samplesPerTrial) / sampleRate;
sampleTimes = [];
for iTrial = 1:originalTrialCount
    sampleTimes = [sampleTimes, (trialSampleTimes + ((trialDuration / ...
        originalTrialCount) * (iTrial - 1)))];
end


end