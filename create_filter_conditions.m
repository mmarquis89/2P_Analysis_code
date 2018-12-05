function [conditionFilters, conditionNames] = create_filter_conditions(primaryFiltName, allFilts, allFiltNames)
%========================================================================================================================= 
% CREATES LISTS OF BEHAVIORAL FILTER COMBINATIONS
%
% Creates lists of combinations of a set of filter conditions for a specificed primary event type.
%
% INPUTS:
%       primaryFiltName = a string containing the name of the primary filter (i.e. event type)
%
%       allFilts = 
%
%       allFiltNames = 
%
% OUTPUTS:
%       conditionFilters = 
%
%       conditionNames = 
%
%
%==========================================================================================================================



% Separate out main filter from secondary filteres
secondaryFilts = allFilts;
mainFiltName = primaryFiltName;
secondaryFiltNames = allFiltNames;

% Get index vector of all possible secondary filter combinations
for iFilt = 1:numel(secondaryFiltNames)
    filtInds{iFilt} = 1:numel(secondaryFiltNames{iFilt});
end
combInds = combvec(filtInds{:})';

% Create condition names and filter arrays
conditionFilters = [];
conditionNames = [];
for iComb = 1:size(combInds, 1)
    currName = '';
    currFilt = [];
    for iFilt = 1:size(combInds, 2)
        currName = [currName, '_', secondaryFiltNames{iFilt}{combInds(iComb, iFilt)}];
        currFilt = [currFilt; secondaryFilts{iFilt}(combInds(iComb, iFilt), :)];
    end
    conditionNames{iComb} = [mainFiltName, currName];
    conditionFilters{iComb} = currFilt;
end
conditionNames = conditionNames';
conditionFilters = conditionFilters';

end