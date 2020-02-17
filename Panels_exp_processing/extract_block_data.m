function bl = extract_block_data(mD, blTrials, varargin)


p = inputParser;
addParameter(p, 'flowSmWin', 30);
addParameter(p, 'MoveThresh', 0.05);
addParameter(p, 'ExpType', 'PB');
addParameter(p, 'NanFlag', 'omitnan');
parse(p, varargin{:});
flowSmWin = p.Results.flowSmWin;
moveThresh = p.Results.MoveThresh;
expType = p.Results.ExpType;

bD = mD(ismember([mD.trialNum], blTrials));

% Check that have compatible values in scalar fields 
if numel(unique([bD.daqSampRate])) > 1 || ...
        numel(unique({bD.expID})) > 1 || ...
        numel(unique([bD.nDaqSamples])) > 1 || ...
        numel(unique([bD.nPlanes])) > 1 || ...
        numel(unique([bD.nPanelsFrames])) > 1 || ...
        numel(unique([bD.nVolumes])) > 1 || ...
        numel(unique([bD.panelsDisplayRate])) > 1 || ...
        numel(unique([bD.trialDuration])) > 1 || ...
        numel(unique([bD.using2P])) > 1 || ...
        numel(unique([bD.volumeRate])) > 1 || ...
        numel(unique([bD.nVolumes])) > 1 || ...
        numel(unique([bD.panelsCycleFrames])) > 1 || ...
        numel(unique([bD.panelsCycleTime])) > 1
    disp('Error: the specified trials are not compatible');
end

% Add data that applies to all trials
bl = [];
bl.daqSampRate = bD(1).daqSampRate;
bl.daqSampTimes = bD(1).daqSampTimes;
bl.expID = bD(1).expID;
bl.nDaqSamples = bD(1).nDaqSamples;
bl.nPanelsFrames = bD(1).nPanelsFrames;
bl.nPlanes = bD(1).nPlanes;
bl.nVolumes = bD(1).nVolumes;
bl.optoStimTiming = bD(1).optoStimTiming;
bl.panelsCycleFrames = bD(1).panelsCycleFrames;
bl.panelsCycleTime = bD(1).panelsCycleTime;
bl.panelsDisplayRate = bD(1).panelsDisplayRate;
bl.panelsFrameTimes = bD(1).panelsFrameTimes;
bl.panelsPattern = bD(1).panelsPattern;
bl.panelsPosX = bD(1).panelsPosX;
bl.panelsPosY = bD(1).panelsPosY;
bl.refImages = bD(1).refImages;
bl.trialDuration = bD(1).trialDuration;
bl.volTimes = bD(1).volTimes;
bl.volumeRate = bD(1).volumeRate;


% Add compiled data
bl.dffArr = []; bl.dffVectAvgRad = []; bl.dffVectAvgWedge = []; bl.dffVectStrength = [];
bl.ftData = []; bl.optoStimOnsetTimes = []; bl.optoStimOffsetTimes = []; bl.rawFlArr = [];
bl.wedgeDffArr = []; bl.wedgeZscoreArr = []; bl.zscoreArr = []; bl.usingOptoStim = [];
bl.trialNum = []; bl.usingPanels = []; bl.wedgeRawFlArr = []; bl.moveDistVols = []; ...
bl.meanVolFlow = []; bl.wedgeExpDffArr = [];
for iTrial = 1:numel(bD)
    bl.dffArr(:, :, iTrial) = bD(iTrial).dffMat;                        % --> [volume, glom, trial]
    bl.rawFlArr(:, :, iTrial) = bD(iTrial).rawFlMat;                    % --> [volume, glom, trial]
    bl.zscoreArr(:, :, iTrial) = bD(iTrial).zscoreMat;                  % --> [volume, glom, trial]
    try % Because some early experiments don't have this one
        bl.expDffArr(:, :, iTrial) = bD(iTrial).expDffMat;              % --> [volume, glom, trial]
    catch; end 
    
    bl.optoStimOnsetTimes{iTrial} = bD(iTrial).optoStimOnsetTimes;      % --> {trial}
    bl.optoStimOffsetTimes{iTrial} = bD(iTrial).optoStimOffsetTimes;    % --> {trial}
    bl.usingOptoStim(iTrial) = bD(iTrial).usingOptoStim;                % --> [trial]
    
    if strcmp(expType, 'PB')
        bl.dffVectAvgRad(:, iTrial) = bD(iTrial).dffVectAvgRad;             % --> [volume, trial]
        bl.dffVectAvgWedge(:, iTrial) = bD(iTrial).dffVectAvgWedge;         % --> [volume, trial]
        bl.dffVectStrength(:, iTrial) = bD(iTrial).dffVectStrength;         % --> [volume, trial]
        bl.wedgeRawFlArr(:, :, iTrial) = bD(iTrial).wedgeRawFlMat;          % --> [volume, wedge, trial]
        bl.wedgeDffArr(:, :, iTrial) = bD(iTrial).wedgeDffMat;              % --> [volume, wedge, trial]
        bl.wedgeZscoreArr(:, :, iTrial) = bD(iTrial).wedgeZscoreMat;        % --> [volume, wedge, trial]
        try
            bl.wedgeExpDffArr(:, :, iTrial) = bD(iTrial).wedgeExpDffMat;    % --> [volume, wedge, trial]
        catch; end
        
    end
    
    bl.trialNum(iTrial) = bD(iTrial).trialNum;                          % --> [trial]
    bl.usingPanels(iTrial) = bD(iTrial).usingPanels;                    % --> [trial]
    
    % Trim FicTrac frames if necessary 
    targetFrames = min(cellfun(@numel, {bD.ftFrameTimes}));
    if numel(bD(iTrial).ftFrameTimes) == targetFrames
       bl.ftFrameTimes = bD(iTrial).ftFrameTimes; 
    end
    % All are --> [frame, trial]
    bl.ftData.intX(:, iTrial) = bD(iTrial).ftData.intX(1:targetFrames);
    bl.ftData.intY(:, iTrial) = bD(iTrial).ftData.intY(1:targetFrames);
    bl.ftData.intHD(:, iTrial) = bD(iTrial).ftData.intHD(1:targetFrames);
    bl.ftData.moveSpeed(:, iTrial) = bD(iTrial).ftData.moveSpeed(1:targetFrames);
    bl.ftData.intFwMove(:, iTrial) = bD(iTrial).ftData.intFwMove(1:targetFrames);
    bl.ftData.intSideMove(:, iTrial) = bD(iTrial).ftData.intSideMove(1:targetFrames);
    bl.ftData.yawSpeed(:, iTrial) = bD(iTrial).ftData.yawSpeed(1:targetFrames);
    bl.ftData.fwSpeed(:, iTrial) = bD(iTrial).ftData.fwSpeed(1:targetFrames);
    bl.ftData.sideSpeed(:, iTrial) = bD(iTrial).ftData.sideSpeed(1:targetFrames);
    
    bl.ftData.moveSpeed(1,iTrial) = bl.ftData.moveSpeed(2,iTrial);
    
    % Use optic flow data to identify movement epochs
    currTrialFlow = bD(iTrial).flowData;
    smFlow = repeat_smooth(currTrialFlow, 20, 'dim', 2, 'smwin', flowSmWin);
    smFlow = smFlow - min(smFlow(:));
    moveFrames = smFlow > moveThresh;
   
    % For each quiescence frame, find the distance to the nearest movement frame
    moveFrameDist = zeros(1, numel(moveFrames));
    moveFrameInds = find(moveFrames);
    for iFrame = 1:numel(moveFrames)
        moveFrameDist(iFrame) = min(abs(moveFrameInds - iFrame));
    end
    
    % Convert to volumes
    volFrames = [];
    volTimes = bD(iTrial).volTimes;
    frameTimes = bD(iTrial).flowFrameTimes;
    for iVol = 1:numel(volTimes)
        [~, volFrames(iVol)] = min(abs(volTimes(iVol) - frameTimes));
    end
    bl.moveDistVols(:, iTrial) = moveFrameDist(volFrames) .* ...
            (numel(volTimes) / numel(frameTimes));
    
        
    % Get mean flow data for each volume    
    frameVols = [];
    for iFrame = 1:numel(frameTimes)
       [~, frameVols(iFrame)] = min(abs(frameTimes(iFrame) - volTimes));
    end
    for iVol = 1:numel(volTimes)
        bl.meanVolFlow(iVol, iTrial) = mean(bD(iTrial).flowData(frameVols == iVol));
    end
    
    
end

% Sort fields alphabetically
bl = orderfields(bl);


end