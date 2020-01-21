

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\20200114-3_38A11_ChR_60D05_7f\ProcessedData';
roiFile = 'Behavior_Vid_ROI_Data.mat';


load(fullfile(parentDir, roiFile), 'roiData'); % --> [y, x]
ftVids = dir(fullfile(parentDir, 'FicTrac_video_trial*.mp4'));

meanFlowMags = {};
for iTrial = 1:numel(ftVids)
   currFile = ftVids(iTrial).name;
   disp(currFile);
   tic
   meanFlowMags{iTrial} = optic_flow_calc(fullfile(parentDir, currFile), 'roiMask', roiData);    
   disp(['Flow calculation completed in ', num2str(toc, '%.1f'), ' sec']); 
end

save(fullfile(parentDir, 'flowMags'), 'meanFlowMags')


%% Load and process optic flow data


figure(1);clf;hold on
for i=13%:numel(meanFlowMags)
    plotData = repeat_smooth(meanFlowMags{i}, 20, 'dim', 2, 'smwin', 6);
    plot(plotData - min(plotData));    
end
ylim([0 0.5])