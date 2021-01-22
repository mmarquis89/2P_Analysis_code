parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\gaplessAcq';

if ~isdir(fullfile(parentDir, 'ftFileBackup'))
    mkdir(fullfile(parentDir, 'ftFileBackup'));
end

% Find all FicTrac data in directory
ftFiles = dir(fullfile(parentDir, '*ficTracData.mat'));

for iFile = 1:numel(ftFiles)
    currFileName = ftFiles(iFile).name;
    disp(currFileName)
    
    % Make backup of file
    copyfile(fullfile(parentDir, currFileName), fullfile(parentDir, 'ftFileBackup', currFileName));
    
    % Load current experiment's data
    load(fullfile(parentDir, currFileName), 'ftData');

    for iTrial = 1:numel(ftData)
        
        % Identify locations of any frames that are filled with 'nan' due to issues with videos
        badFrames = isnan(ftData(iTrial).intX);
        
        % Extract data, dropping any nan frames
        intX = ftData(iTrial).intX(~badFrames);
        intY = ftData(iTrial).intY(~badFrames);
        intHD = ftData(iTrial).intHD(~badFrames);
        intFwMove = ftData(iTrial).intFwMove(~badFrames);
        intSideMove = ftData(iTrial).intSideMove(~badFrames);
        moveSpeed = ftData(iTrial).moveSpeed(~badFrames);
        yawSpeed = ftData(iTrial).yawSpeed(~badFrames);
        fwSpeed = ftData(iTrial).fwSpeed(~badFrames);
        sideSpeed = ftData(iTrial).sideSpeed(~badFrames);        
        
        % Use integrated movement data to identify boundaries of original trials
        trialStartFrames = find(intX == 0 & intY == 0 & intHD == 0);
        trialEndFrames = [trialStartFrames(2:end) - 1; numel(intX)];

        % Unwrap heading direction data 
        uwIntHD = unwrap(intHD);

        % Remove resets of all integrated positions to 0 for all fields
        nShortTrials = numel(trialStartFrames);
        for iShort = 2:nShortTrials
            currTrialInds = trialStartFrames(iShort):trialEndFrames(iShort);
            intX(currTrialInds) = intX(currTrialInds) - (intX(trialStartFrames(iShort)) - ...
                    intX(trialEndFrames(iShort - 1)));
            intY(currTrialInds) = intY(currTrialInds) - (intY(trialStartFrames(iShort)) - ...
                    intY(trialEndFrames(iShort - 1)));
            uwIntHD(currTrialInds) = uwIntHD(currTrialInds) - (uwIntHD(trialStartFrames(iShort)) - ...
                    uwIntHD(trialEndFrames(iShort - 1))); 
            intFwMove(currTrialInds) = intFwMove(currTrialInds) - (intFwMove(trialStartFrames(iShort)) - ...
                    intFwMove(trialEndFrames(iShort - 1)));
            intSideMove(currTrialInds) = intSideMove(currTrialInds) - ...
                    (intSideMove(trialStartFrames(iShort)) - intSideMove(trialEndFrames(iShort - 1)));
        end

        % Re-wrap the corrected heading direction data
        intHD = mod(uwIntHD, (2*pi));

        % Correct the speeds by averaging across gaps of original trials
        moveSpeed(trialStartFrames(2:end)) = (moveSpeed(trialStartFrames(2:end) - 1) + ...
                moveSpeed(trialStartFrames(2:end) + 1)) ./ 2;     
        yawSpeed(trialStartFrames(2:end)) = (yawSpeed(trialStartFrames(2:end) - 1) + ...
                yawSpeed(trialStartFrames(2:end) + 1)) ./ 2; 
        fwSpeed(trialStartFrames(2:end)) = (fwSpeed(trialStartFrames(2:end) - 1) + ...
                fwSpeed(trialStartFrames(2:end) + 1)) ./ 2; 
        sideSpeed(trialStartFrames(2:end)) = (sideSpeed(trialStartFrames(2:end) - 1) + ...
                sideSpeed(trialStartFrames(2:end) + 1)) ./ 2; 
            
        % Replace the orginal vectors with the corrected ones, putting the 'nan' frames back 
        ftData(iTrial).intX = nan; 
        ftData(iTrial).intY = nan;
        ftData(iTrial).intHD = nan;
        ftData(iTrial).intFwMove = nan;
        ftData(iTrial).intSideMove = nan;
        ftData(iTrial).moveSpeed = nan;
        ftData(iTrial).yawSpeed = nan;
        ftData(iTrial).fwSpeed = nan;
        ftData(iTrial).sideSpeed = nan;        
        ftData(iTrial).intX(~badFrames) = intX;
        ftData(iTrial).intY(~badFrames) = intY;
        ftData(iTrial).intHD(~badFrames) = intHD;
        ftData(iTrial).intFwMove(~badFrames) = intFwMove;
        ftData(iTrial).intSideMove(~badFrames) = intSideMove;
        ftData(iTrial).moveSpeed(~badFrames) = moveSpeed;
        ftData(iTrial).yawSpeed(~badFrames) = yawSpeed;
        ftData(iTrial).fwSpeed(~badFrames) = fwSpeed;
        ftData(iTrial).sideSpeed(~badFrames) = sideSpeed;

    end%iTrial
    
    % Save corrected data for this experiment
    save(fullfile(parentDir, currFileName), 'ftData');
    
    

end
