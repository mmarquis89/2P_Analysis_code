
function cycleData = get_bar_cycle_data(trialData, flSmWin)
%===================================================================================================
% For a trial with an open loop swinging bar (or other shape) visual stimulus, generates a table 
% containing timing, visual stim position, and fluorescence data for each individual cycle of the 
% visual stim. Also calculates the cycle's vector strength and phase for each wedge or glomerulus 
% ROI. 
% 
% Note: incomplete cycles that were cut off at the end of a trial are included, but they are marked 
% as incomplete in the "fullCycles" column and their vector strength/phase are filled with NaNs.
%
% INPUTS:
%
%       trialData   = source data and metadata table for a single trial. Will use the following 
%                     specific fields: 
%                               nPanelsFrames
%                               panelsDisplayRate
%                               nVolumes
%                               volumeRate
%                               panelsPosX
%                               rawFl
%                               trialDff
%                               expDff
% 
%       flSmWin     = width of the gaussian smoothing window to use on the fluorescence data before 
%                     calculating each cycle's vector strength and phase.
% 
% OUTPUTS:
%
%       cycleData   = table containing data for each individual bar cycle. The first three fields 
%                     uniquely identify each row. The field names are:
%                          expID               
%                          trialNum
%                          cycleNum            (cycle number relative to the start of the trial)
%                          cycStartVol         (trial volume when the current cycle starts)
%                          cycEndVol           (trial volume when the current cycle ends)
%                          fullCycle           (whether cycle was full or cut off at end of trial)
%                          cycBarPos           (index of the bar position during each volume)
%                          cycBarPhase         (same as above, but scaled from 0-2*pi)
%                          cycVolTimes         (the times of each volume in the cycle)
%                          trialVolTimes       (times of each volume relative to start of trial)
%                          cycRawFl            (matrix of raw fluorescence for each volume and ROI)
%                          cycTrialDff         (same as above, but for dF/F with a trial baseline)
%                          cycExpDff           (same as above, but for dF/F with a exp. baseline)
%                          cycVectorStrength   (current cycle's vector strength value for each ROI)
%                         cycVectorPhase      (current cycles vector phase value for each ROI)
% 
%===================================================================================================


td = trialData;

% Identify the panels X location during each volume
panelsFrameTimes = double(1:td.nPanelsFrames) ./ td.panelsDisplayRate;
volTimes = (1:td.nVolumes) ./ td.volumeRate;
panelsPosVols = [];
for iVol = 1:td.nVolumes
    [~, currVol] = min(abs(panelsFrameTimes - volTimes(iVol)));
    panelsPosVols(iVol) = td.panelsPosX{:}(currVol);
end

% Convert panels position to a phase
panelsBarPhase = 2*pi * (panelsPosVols ./ max(panelsPosVols));

% Identify the start and end indices of each bar cycle
panelsPosStr = regexprep(num2str(panelsPosVols == 0), ' ', '');
cycStartVols = [1, regexp(panelsPosStr, '(?<=0)1')];
cycEndVols = [cycStartVols(2:end) - 1, td.nVolumes];

% Drop last cycle of the trial if it is incomplete
fullCycles = ones(size(cycStartVols))';
if panelsPosVols(cycEndVols(end)) ~= panelsPosVols(cycEndVols(1))
    fullCycles(end) = 0;
%     cycStartVols = cycStartVols(1:end-1);
%     cycEndVols = cycEndVols(1:end-1);
end
nCycles = numel(cycStartVols);

% Create a table with data for each bar cycle
cycleData = [];
for iCycle = 1:nCycles
    
    % Create new row
    newRow = table(td.expID, td.trialNum, iCycle, 'VariableNames', {'expID', 'trialNum', ...
            'cycleNum'});
    
    % Identify current cycle's start and end volumes
    sVol = cycStartVols(iCycle);
    eVol = cycEndVols(iCycle);
        
    % Add bar position and timing info
    newRow.cycStartVol = sVol;
    newRow.cycEndVol = eVol;
    newRow.fullCycle = fullCycles(iCycle);
    newRow.cycBarPos = {panelsPosVols(sVol:eVol)'};
    newRow.cycBarPhase = {panelsBarPhase(sVol:eVol)'};
    newRow.cycVolTimes = {(1:numel(newRow.cycBarPos{:}))' ./ td.volumeRate};
    newRow.trialVolTimes = {volTimes(sVol:eVol)'};
    
    % Add fluorescence data for the current cycle
    newRow.cycRawFl = {td.rawFl{:}(sVol:eVol, :)};
    newRow.cycTrialDff = {td.trialDff{:}(sVol:eVol, :)};
    newRow.cycExpDff = {td.expDff{:}(sVol:eVol, :)};
    
    % Use expDff to calculate vector strength and phase for each wedge or glomerulus
    smDff = smoothdata(newRow.cycExpDff{:}, 1, 'gaussian', flSmWin);
    vStrength = [];
    vPhase = [];
    for iRoi = 1:size(smDff, 2)
        [vStrength(iRoi), vPhase(iRoi)] = calculate_vector_strength(newRow.cycBarPhase{:}, ...
                smDff(:, iRoi));
    end
    if fullCycles(iCycle)
        newRow.cycVectorStrength = {vStrength};
        newRow.cycVectorPhase = {vPhase};
    else
        newRow.cycVectorStrength = {nan(size(vStrength))};
        newRow.cycVectorPhase = {nan(size(vPhase))};
    end
    
    
    cycleData = [cycleData; newRow];
end%iCycle

end%function











