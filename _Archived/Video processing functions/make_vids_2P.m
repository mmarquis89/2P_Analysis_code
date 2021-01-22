function msg = make_vids_2P(sid, parentDir, frameRate)
%==============================================================================================================================
% CREATE MOVIES FROM .TIF FILES FROM 2P EXPERIMENT
% Creates an .avi movie for each trial of a 2P imaging experiment from the .tif files captured by the fly behavior camera.
% The videos are saved in a folder within the parent directory named '_Movies'.
%
% Inputs:
%       sid       = the session ID of the videos you want to process.
%       parentDir = the directory containing the subfolders with videos from all sessions.
%                   e.g. 'U:\2P Behavior Video\2017_07_30' with subfolder e.g. '\bdata_OdorCenterWind_20170730_170018_sid_0_tid_0'
%       frameRate = the frame rate that the behavior camera was acquiring at.
%===============================================================================================================================

savePath = fullfile(parentDir, '_Movies');
dirFiles = dir(fullfile(parentDir, ['*sid_', num2str(sid), '_t*']));
sessionDirs = dirFiles([dirFiles.isdir]); % Drop any non-folder items in the directory (e.g. ball/metadata files)

disp('Creating videos...')

    frameCounts = [];
    for iTrial = 1:length(sessionDirs)
        
        % Pad the trial number with leading zeros if necessary to ensure correct filename sorting
        if iTrial < 10
            padStr = '00';
        elseif iTrial < 100
            padStr = '0';
        else
            padStr = '';
        end
        
        % Initialize variables
        trialStr = ['sid_', num2str(sid), '_tid_', padStr, num2str(iTrial)];
        disp(trialStr)
        folderName = sessionDirs(iTrial).name;
        sourcePath = fullfile(parentDir, folderName);
        
        % Create save directory if it doesn't already exist
        if ~isdir(savePath)
            mkdir(savePath);
        end
        
        % Make sure this video doesn't already exist
        if exist(fullfile(savePath, [trialStr, '.mp4']), 'file') == 0
            
            % Create video writer object using MPEG-4 H.264 compression for Anvil compatibility
            outputVid = VideoWriter(fullfile(savePath, [trialStr, '.mp4']), 'MPEG-4');
            outputVid.FrameRate = frameRate;
            open(outputVid)
            
            % Make sure there's at least one image file in this trial's directory
            currFiles = dir(fullfile(parentDir, folderName, '*.tif'));
            if ~isempty(currFiles)
                currFrames = sort({currFiles.name}');
                
                % Write each .tif file to video
                for iFrame = 1:length(currFrames)
                    
                    % Read image
                    currImg = imread(fullfile(sourcePath, currFrames{iFrame}));
                    
                    % Write frame to video
                    writeVideo(outputVid, currImg);                    
                end
                
                % Record number of frames in log
                frameCounts(iTrial).nFrames = length(currFrames);
                frameCounts(iTrial).trial = trialStr;
            else
                % Record number of frames in log
                frameCounts(iTrial).nFrames = 0;
                frameCounts(iTrial).trial = trialStr;
            end%if
            
            close(outputVid)
        else
            disp(['Video already exists...skipping ', trialStr]);
        end%if
    end%iTrial
    
    % Save frame count log
    save(fullfile(savePath, ['sid_', num2str(sid) '_frameCountLog.mat']), 'frameCounts')    
    msg = 'Videos created successfully!';

end%function