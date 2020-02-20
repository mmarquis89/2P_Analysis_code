function sampleInds = time2idx(queryTimes, referenceTimes)

    sampleInds = [];
    for iTime = 1:numel(queryTimes)
        [~, sampleInds(iTime)] = min(abs(referenceTimes - queryTimes(iTime)));
    end
   

end