outputTbl = [];
for iType = 1:size(neuronList, 1)
    currData = partnerTypeTbl(strcmp(partnerTypeTbl.neuronID, neuronList.neuronID{iType}), :);
    currData = sortrows(currData, 'synapseCount', 'descend');
    
    nPartners = size(currData, 1);
    pctCutoffInd = round(nPartners * .1);
%     pctCutoffInd = 10;

    topPartnerTotalPct = cumsum(currData.connectionWeightPct(1:pctCutoffInd));
    cumPctVal = topPartnerTotalPct(end);
    cutoffPartnerVal = currData.connectionWeightPct(pctCutoffInd);
    
    topPartnerMeanPct = mean(currData.connectionWeightPct(1:pctCutoffInd));
    
    
    
    outputTbl = [outputTbl; table(neuronList.neuronID(iType), neuronList.neuronName(iType), ...
            cumPctVal, topPartnerMeanPct, cutoffPartnerVal, 'variableNames', {'neuronID', ...
            'neuronName', 'topPartnerCumPct', 'topPartnerMeanPct', 'cutoffPercentile'})];
        
        
end
outputTbl = sortrows(outputTbl, 5);