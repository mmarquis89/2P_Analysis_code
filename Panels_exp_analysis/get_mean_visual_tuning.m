function outputTbl = get_mean_visual_tuning(tbl, timeWindows)
%===================================================================================================
% Accepts a table (timeWindows) specifying one or more time windows in which to extract mean visual 
% tuning data from a source data table with one row for each trial, and returns the results in an 
% output table with one row for each time window that was requested.
%
% The visual tuning data is rotated so that the center position is directly in front of the fly (the
% "barAngles" field gives the exact angles relative to the fly for each of the 96 positions).
%
% 
% INPUTS:
% 
%       tbl             = the source data/metadata table to extract the tuning data from. Needs to 
%                         have a row for each expID+trialNum combination in the timeWindows table.
%                         the specific fields that will be used are:
%                               expDff
%                               trialNum
%                               volTimes
%                               nPanelsFrames
%                               panelsDisplayRate
%                               panelsPosX
%                               nVolumes
%                               rawFl
%                               trialDff
%                               expDff
% 
%       timeWindows     = a table specifying the time windows to get tuning data from, with columns:
%                               [expID][trialNum][startTime][endTime]
%                         Start and end times should be given in seconds.  
% 
% OUTPUTS:
% 
%       outputTbl       = a copy of the timeWindows input table, but with extra columns containing 
%                         the list of bar angles relative to the fly (in deg) and the mean fl data 
%                         for each position. Output column names are: 
%                               [barAngles][rawFlTuning][trialDffTuning][expDffTuning]
%                           
%
%===================================================================================================

outputTbl = [];
for iWin = 1:size(timeWindows, 1)
    
    % Get correct source data row and time window timing info
    currTbl = tbl(strcmp(tbl.expID, timeWindows.expID{iWin}) & ...
            tbl.trialNum == timeWindows.trialNum(iWin), :);
    startTime = timeWindows.startTime(iWin);
    endTime = timeWindows.endTime(iWin);
    volTimes = currTbl.volTimes{:};
    
    % Identify the panels X location during each volume
    panelsFrameTimes = double(1:currTbl.nPanelsFrames) ./ currTbl.panelsDisplayRate;
    panelsPosX = currTbl.panelsPosX{:};
    panelsPosVols = [];
    for iVol = 1:currTbl.nVolumes
        [~, currVol] = min(abs(panelsFrameTimes - volTimes(iVol)));
        panelsPosVols(iVol) = panelsPosX(currVol);
    end
    
    % Extract fluorescence and bar location data for the current time window
    currWinVols = volTimes > startTime & volTimes < endTime;
    currWinPanelsPos = panelsPosVols(currWinVols);
    rawFl = currTbl.rawFl{:}(currWinVols, :);
    trialDff = currTbl.trialDff{:}(currWinVols, :);
    expDff = currTbl.expDff{:}(currWinVols, :);
    
    % Get mean fluorescence at each bar location
    barPositions = (0:95)';
    rawFlTuning = []; trialDffTuning = []; expDffTuning = [];
    for iPos = 1:numel(barPositions)
        currPos = barPositions(iPos);
        rawFlTuning(iPos, :) = mean(rawFl(currWinPanelsPos == currPos, :), 1, 'omitnan');
        trialDffTuning(iPos, :) = mean(trialDff(currWinPanelsPos == currPos, :), 1, 'omitnan');
        expDffTuning(iPos, :) = mean(expDff(currWinPanelsPos == currPos, :), 1, 'omitnan');
    end
    
    % Shift data so the center is directly in front of the fly
    rawFlTuning = rawFlTuning([93:96, 1:92], :);
    trialDffTuning = trialDffTuning([93:96, 1:92], :);
    expDffTuning = expDffTuning([93:96, 1:92], :);
    barAngles = (-180 + 3.75):3.75:180;
    
    % Append tuning data to output table
    newRow = timeWindows(iWin, :);
    newRow.barAngles = {barAngles};
    newRow.rawFlTuning = {rawFlTuning};
    newRow.trialDffTuning = {trialDffTuning};
    newRow.expDffTuning = {expDffTuning};
    outputTbl = [outputTbl; newRow];    
    
end%iWin

end%function
