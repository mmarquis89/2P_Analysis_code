
startDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';
parentDir = uigetdir(startDir, 'Select an experiment directory');
roiFile = 'Behavior_Vid_ROI_Data.mat';


load(fullfile(parentDir, roiFile), 'roiData'); % --> [y, x]
% ftVids = dir(fullfile(parentDir, 'FicTrac_video_trial*.mp4'));
ftVids = dir(fullfile(parentDir, '*.avi'));

meanFlowMags = {};
for iTrial = 1:numel(ftVids)
   currFile = ftVids(iTrial).name;
   disp(currFile);
   tic
   meanFlowMags{iTrial} = optic_flow_calc(fullfile(parentDir, currFile), 'roiMask', roiData);    
   disp(['Flow calculation completed in ', num2str(toc, '%.1f'), ' sec']); 
end

save(fullfile(parentDir, 'flowMags'), 'meanFlowMags')
