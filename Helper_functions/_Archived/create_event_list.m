function eventList = create_event_list(onsetArr, offsetArr)
%===================================================================================================
% 
% Creates a list of each event in an array in 2D annotation format (see the documentation for 
% annot2lin and lin2annot for more details). Each event must have both an onset and an offset. The 
% starting and ending indices of each event will be listed vertically in the output eventList.
%
% INPUTS:
%       onsetArr  = logical array in 2D annotation format (see annot2lin description) containing the
%                   locations of all event onsets.
%
%       offsetArr = same thing, but for the event offsets
%
% OUTPUTS:
%       eventList = nEvents x 3 array with columns: [onsetIdx, offsetIdx, trialNum]
%
%===================================================================================================

eventList = [];
for iTrial = 1:size(onsetArr, 1)
        currTrialOnsets = find(onsetArr(iTrial,:));
        currTrialOffsets = find(offsetArr(iTrial,:));
        for iEvent = 1:length(currTrialOnsets)
            eventList(end+1, :) = [currTrialOnsets(iEvent), currTrialOffsets(iEvent), iTrial];
        end
end

end