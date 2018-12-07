%==================================================================================================
%%% EXTRACT SAMPLE FRAMES FROM ANATOMY STACKS -----------------------------------------------------
%% ==================================================================================================

expDates = {...
    '2018_04_26_exp_1'
            }

targetPlanes = [ 200 ...
    ];

fileStr = '*Stack_*.tif';

for iExp = 1:length(expDates)
    
    expDate = expDates{iExp}
    targetPlane = targetPlanes(iExp);
    
    dirPath = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate];
    tic
    disp(['Extracting sample frames from anatomy stacks...']);
    extract_sample_frames(dirPath, fileStr, targetPlane);
    writeToLog(sprintf('%s sample frames extracted in %s min', expDate, num2str(round(toc/60, 1))));
    disp(['Extracting sample frames took ', num2str(round(toc/60, 1)) ' min']);
end
clear expDates targetPlanes fileStr dirPath
% -------------------------------------------------------------------------------------------------

%% ===================================================================================================
%%% SELECT ROI FOR OPTIC FLOW CALCULATION
%% ===================================================================================================

expDates = {...
    '2018_11_11_exp_3' ...
            }
     
 for iExp = 1:length(expDates)
     
     expDate = expDates{iExp}
     
     %%% Define ROIs for optic flow combined vids
     parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate];
     select_video_ROIs(parentDir);
 end
 
%% ==================================================================================================
%%% PROCESS ANVIL ANNOTATION DATA ------------------------------------------------------------------
%% ==================================================================================================

FRAME_RATE = 25;

trialDuration = 300;
expDate = '2018_10_30_exp_1';
sid = 0;

parentDir = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, '\_Movies'];
saveDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)];
annotationFileName = [expDate, '_Annotation.txt'];
% annotationFileName = [expDate, '_sid_', num2str(sid), '_Annotation.txt'];

tic
process_anvil_annotations(sid, parentDir, saveDir, annotationFileName, FRAME_RATE, trialDuration);
disp(['Processing Anvil annotations took ', num2str(round(toc/60, 1)) ' min']);

clear parentDir saveDir annotationFileName
% -------------------------------------------------------------------------------------------------


%% ==================================================================================================
%%% ARCHIVE FILES-----------------------------------------------------------------------------------
%% ==================================================================================================

% -------------------------------------------------------------------------------------------------
% To-do:        
%        

expDateTemp = '2018_10_09_exp_1';
parentDirTemp = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDateTemp];
sidTemp = 0;

%%% Archive raw anatomy stacks
archiveName = 'AnatomyStacks';
filterString = '*Stack*';
system7zip(parentDirTemp, archiveName, '7z', filterString, 1)

%%% Archive raw imaging data files
archiveName = ['TrialData_sid_', num2str(sidTemp)];
filterString = ['*sid_', num2str(sidTemp), '_b*'];
system7zip(parentDirTemp, archiveName, '7z', filterString, 1)

%%% Archive raw video frames
parentDirTemp = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDateTemp];
archiveName = ['sid_', num2str(sidTemp), '_RawFrames'];
system7zip(parentDirTemp, archiveName, '7z', ['*sid_', num2str(sidTemp), '_b*'], 1);

clear parentDir archiveName filterString
% -------------------------------------------------------------------------------------------------