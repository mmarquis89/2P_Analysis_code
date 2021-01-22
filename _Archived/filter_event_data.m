function [filtOutput] = filter_event_data(eventList, filterEventData, analysisWindowIn, filterEventWindows, filterDirections, ...
                                        varargin)
%========================================================================================================================
% 
% This function takes a list of annotation events (e.g. odor delivery, ball stopping, locomotion) and 
% determines whether they meet a list of criteria based on whether other types of events ("filter" events)
% occurred just before, during, or after the target event. The number of filter event types, size of the 
% evaluation time windows, and direction of the filtering (i.e. whether the filter event is required to be 
% present or absent to match the filter) are all adjustable.
%
% INPUTS:
%       eventList           = an M x 3 array specifying the onset volume, offset volume, and trial number of 
%                             an event (as columns, in that order)
%
%       filterEventData     = an A x B x N array containing the annotation data (in volumes) for each type of 
%                             event that you want to use to filter the target events in eventList. The 
%                             dimensions of the array should be: [trial, volume, eventType]
%
%       analyisWindowIn     = a 1 x 2 array specifying the length (in seconds) of the [baseline, response] 
%                             periods you will use in your analysis of the target event. The filter windows will 
%                             be applied starting at either side of this window.
% 
%       filterEventWindows  = an N x 2 array specifying the [pre-stim, post-stim] time windows in seconds 
%                             around the analysis window to use for filtering. Should have one row for each 
%                             type of filter event.
%
%       filterDirections    = an N x 3 array of numeric values indicating whether each filter type should be 
%                             inclusive or exclusive during the pre-event, during response, and post-response
%                             filter windows. A value of 1 means inclusive, -1 means exclusive, and 0 means that you
%                             don't care whether the filter event happened during that period. For example, 
%                             passing [-1 -1 -1 ; -1 1 0] would match the target events for which:
%                                - The first type of filter event didn't occur at any point within the filter window.
%                                - The second type of filter event occurred during the target event's response 
%                                  period, but not in the window before the event (the post-response window is 
%                                  irrelevant in this case because the filter direction is zero).
%
%       'volumeRate'          = (default = 6.44) imaging volume acquisition rate, for converting from seconds to volumes
%
%       'offsetAlign'         = <OPTIONAL> boolean value specifying whether or not to align the filtering to the onset or 
%                             the offset of each target event (default == 0).
%
%       'overshoot'           = <OPTIONAL> a boolean value specifying whether to allow events that are shorter in duration 
%                             than the response (if offsetAlign == 0) or the baseline (if offsetAlign == 1) period. 
%                             Defaults to zero if no value is provided. 
%
% OUTPUTS:
%       filtOutput          = a logical vector indicating whether each event in the eventList passed the filtering 
%                             process for all the filter event types.
% 
%       skipCount           = the number of events that were discarded due to overshoot or proximity to trial start 
%                             or end
%
%=====================================================================================================================
%
% EXAMPLE of the relative timing of different components of the filtering when aligned to a target event onset:
%
%       filtWin (AKA fW)        = [2, 1]
%       analysisWin (AKA aW)    = [1, 1]
%       filtDirections (AKA fD) = [-1  1  1]
%       alignVol                = startVol
%
%                                                 |------------event------------|          
%       <--------fW(1)--------><----aW(1)---->[alignVol]<----aW(2)----><----fW(2)---->
%       [------------------fD(1)-----------------][-------fD(2)-------][----fD(3)----]
% 
% Based on the filter directions provided here the filter would identify all target events that are 1) preceded by 
% at least 3 seconds with no instances of the filter event, 2) have the filter event occurring at some 
% point within the response period, and 3) have the filter event continuing throughout the end of the 
% filter window.
%
%
% Here's the same example if offsetAlign == 1:
%
%                 |-------------event-------------|           
%       <--------fW(1)--------><----aW(1)---->[alignVol]<----aW(2)----><----fW(2)---->
%       [------------------fD(1)-----------------][-------fD(2)-------][----fD(3)----]
%
% 
% If overshoot == 1, the following target event would not pass the filtering because the event ends before 
% making it to the end of the analysis window (this would also fail if offsetAlign == 1 because it would not
% reach the beginning of the analysis window: 
% 
%                                                 |----event----|           
%       <--------fW(1)--------><----aW(1)---->[alignVol]<----aW(2)----><----fW(2)---->
%       [------------------fD(1)-----------------][-------fD(2)-------][----fD(3)----]
%
%==============================================================================================================

% Parse optional arguments
p = inputParser;
addOptional(p, 'optArg', []);
addParameter(p, 'volumeRate', 6.44);
addParameter(p, 'offsetAlign', 0);
addParameter(p, 'overshoot', 0);
parse(p, varargin{:});
volumeRate = p.Results.volumeRate;
offsetAlign = p.Results.offsetAlign;
overshoot = p.Results.overshoot;
if ~isempty(p.Results.optArg)
   if p.Results.optArg == 1
      offsetAlign = p.Results.optArg; % If last argument has no name/value pair, still assume it is offsetAlign if == 1
   end
end

% Calculate constants
nVolumes = size(filterEventData, 2);
nEvents = size(eventList, 1);
nFilters = size(filterEventData, 3);

% Convert from seconds to volumes
analysisWindow = floor(analysisWindowIn .* volumeRate);
filterEventWindows = floor(filterEventWindows .* volumeRate);

% Fill in a logical array for each event and filter
filtArr = zeros(nEvents, nFilters);
skipArr = filtArr;
for iFilt = 1:nFilters
    filtWin = filterEventWindows(iFilt, :);
    for iEvent = 1:nEvents
        
        % Get target event information
        [eventStart, eventEnd, trial] = split_vector(eventList(iEvent, :));
        if offsetAlign
            alignVol = eventEnd;
        else
            alignVol = eventStart;
        end
  
        % Calculate the starting and ending volumes of the filtering window
        startVol = alignVol - analysisWindow(1) - filtWin(1);
        endVol = alignVol + analysisWindow(2) + filtWin(2);
        
        % Skip any events that are shorter than the analysis period unless overshoot == 1
        eventDur = eventEnd - eventStart;
        if offsetAlign
            overshootVols = analysisWindow(1);
        else
            overshootVols = analysisWindow(2);
        end
        if ~overshoot && (eventDur < overshootVols)
            skipArr(iEvent, :) = 1;
            continue
        end
        
        % Force non-match if the start or end of the trial is within the analysis filter windows
        if ((alignVol - analysisWindow(1)) < 1) || ((alignVol + analysisWindow(2)) >= nVolumes)
            skipArr(iEvent, :) = 1;
            continue
        end
               
        % Change start/end vols if the beginning or end of the trial falls in the filter windows
        if (startVol < 1) && ((startVol + filtWin(1)) > 0)
            startVol = 1;
        end
        if (endVol > nVolumes) && ((endVol - filtWin(2)) < nVolumes)
            endVol = nVolumes;
        end
        
        % Determine whether any filter events occurred within each filtering period
        if sum(filterEventData(trial, startVol:(alignVol - 1), iFilt))
            preFilt = 1;
        else
            preFilt = -1;
        end
        if sum(filterEventData(trial, alignVol:(alignVol + analysisWindow(2) - 1), iFilt))
            stimFilt = 1;
        else
            stimFilt = -1;
        end
        if sum(filterEventData(trial, ((alignVol + analysisWindow(2)):endVol), iFilt))
            postFilt = 1;
        else
            postFilt = -1;
        end
        
        % Check whether the current event matches the filter
        filterVec = [preFilt, stimFilt, postFilt];
        filterVec(filterDirections(iFilt, :) == 0) = 0; % Force matching for windows you don't care about
        if sum(filterVec == filterDirections(iFilt, :)) == 3
            filtArr(iEvent, iFilt) = 1;
        end

    end%for iEvent
    
end%for iFilt


% Force non-matches for any skipped trials 
filtArr(logical(skipArr)) = 0;
skipCount = sum(skipArr(:,1));
disp(['Skipped ' num2str(skipCount), ' of ', num2str(nEvents), ' events due to overshoot settings or proximity to trial start/end'])

filtOutput = logical(sum(filtArr, 2) == nFilters);

end