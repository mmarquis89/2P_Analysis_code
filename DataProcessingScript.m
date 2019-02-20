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
    '2018_12_18_exp_2' ...
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

expDateTemp = '2018_10_30_exp_1';
parentDirTemp = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDateTemp];
sidTemp = 0;

%%% Archive raw anatomy stacks
parentDirTemp = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDateTemp];
archiveName = 'AnatomyStacks';
filterString = '*Stack*';
system7zip(parentDirTemp, archiveName, '7z', filterString, 1)

%%% Archive raw imaging data files
parentDirTemp = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDateTemp];
archiveName = ['TrialData_sid_', num2str(sidTemp)];
filterString = ['*sid_', num2str(sidTemp), '_b*'];
system7zip(parentDirTemp, archiveName, '7z', filterString, 1)

%%% Archive raw videos
parentDirTemp = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDateTemp];
archiveName = ['sid_', num2str(sidTemp), '_RawVids'];
system7zip(parentDirTemp, archiveName, '7z', ['*sid_', num2str(sidTemp), '_b*'], 1);

%%% Archive block vid frames
parentDirTemp = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDateTemp];
archiveName = ['sid_', num2str(sidTemp), '_blockVids'];
system7zip(parentDirTemp, archiveName, '7z', 'BlockVids', 1);

%%% Archive individual vid frames
parentDirTemp = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDateTemp];
archiveName = ['sid_', num2str(sidTemp), '_RawVids'];
system7zip(parentDirTemp, archiveName, '7z', '20181030*', 1);

%% Archive all 2017 behavior vid files


% parentDirTemp = ['D:\Dropbox (HMS)\2P Data\Behavior Vids\2017'];
% allExpDirs = dir(fullfile(parentDirTemp, '2017*'));

% dirNames = string({allExpDirs.name});


% system7zip(parentDirTemp, dirNames{1}, '7z', dirNames{1}, 1);
% dirNames(1) = [];
















%% -------------------------------------------------------------------------------------------------

parentDirTemp = ['D:\Dropbox (HMS)\2P Data\Imaging Data\2017'];
dirName = '2017_10_26_TH_Gal4_exp_2';
expDir = fullfile(parentDirTemp, dirName);

% Create folder for compiled analysis data
if ~isdir(fullfile(parentDirTemp, 'analysis_dirs'))
   mkdir(fullfile(parentDirTemp, 'analysis_dirs')) 
end

% Get names of any sid dirs 
sidDirs = dir(fullfile(expDir, '*sid_*'));
sidDirs = sidDirs([sidDirs.isdir]);

for iSid = 1:numel(sidDirs)
   
    % Check whether the sid dir contains an analysis folder
    analysisDir = dir(fullfile(expDir, sidDirs(iSid).name, 'a*lysis'));
    if ~isempty(analysisDir)
       % Copy analysis data before archiving experiment
       copyfile(fullfile(analysisDir.folder, analysisDir.name), fullfile(parentDirTemp, 'analysis_dirs', dirName, sidDirs(iSid).name));
    end
end

% Archive entire exp dir
system7zip(parentDirTemp, dirName, '7z', dirName, 1);


