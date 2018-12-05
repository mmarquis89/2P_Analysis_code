function [outputMetadata, wholeSession] = load_imaging_data(parentDir, sessionDataFile, varargin)
%=======================================================================================================================================================
%
%  Loads a data file containing 2P imaging data. Also does some basic
%  pre-processing/metadata extraction before returning it all as a single data structure.
%
%
%       parentDir                = the path to a directory containing all the data/metadata files
%
%       sessionDataFile          = a session data .mat file containing just an array of imaging data
%                                  named 'wholeSession' with dimensions [y, x, plane, volume, trial]
%
%  OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%
%  These are collectively the default assumtions for filenames, which should all be in the same directory as the
%  selected session data file. Pass arguments to override, or [] to prompt a file selection dialog for the
%  indicated filename.
%
%       'sid'                    = (default: looks for it in sessionDataFile name)
%       'AnnotFile'              = (default: 'Annotations.mat')
%       'imgMetadataFile'        = (default: 'imgMetadata.mat')
%       'refImgFile'             = (default: 'refImages_Reg.mat')
%       'LoadSessionData'        = (default: 0) Specifies whether or not to return the entire session array
%
%  OUTPUT:
%       outputData = a structure with the following fields:
%               behaviorLabels   = the labels corresponding to each number in the behavior annotation data
%               expDate          = date of experiment in YYYY_MM_DD format
%               goodTrials       = logical vector specifying all the trials that are not missing any behavior video frames
%               laserPower       = power of the 2P laser during imaging
%               nFrames          = the number of frames/trial in the behavior video/annotation data
%               nPlanes          = the total number of imaging planes in the session
%               nTrials          = the total number of trials in the session
%               nVolumes         = the number of imaging volumes per trial in the session
%               origFileNames    = cell array containing the file names of all the raw .tif files that generated the input data structure
%               refImg           = 1 x nPlanes cell array, with each cell containing the average of one plane across all volumes and trials [y, x]
%               scanFrameRate    = the acquisition frame rate of the 2P data
%               sid              = the session ID of the data that was loaded
%               stimDurs         = 1 x nTrials numeric vector containing the duration of the stimulus in each trial
%               stimOnsetTimes   = 1 x nTrials numeric vector containing the onset time of the stimulus in seconds for each trial
%               stimSepTrials    = a structure with 1 x nTrials logical vectors for each stimType
%               stimTypes        = cell array containing the names of each unique trialType in the session
%               trialAnnotations = cell array with behavior annotation data (see process_anvil_annotations() for more info)
%               trialDuration    = total duration of trial in seconds
%               trialType        = cell array with the stimulus type for each trial that went in the data structure
%               volumeRate       = volume acquisition rate for the imaging data
%
%      wholeSession =  the imaging data array with dimensions [y, x, plane, volume, trial].
%
%========================================================================================================================================================
try
% Parse optional arguments
p = inputParser;
addParameter(p, 'sid', []);
addParameter(p, 'StimMetadataFile', 'stimMetadata.mat');
addParameter(p, 'AnnotFile', 'Annotations.mat');
addParameter(p, 'imgMetadataFile', 'imgMetadata.mat');
addParameter(p, 'refImgFile', 'refImages_Reg.mat')
addParameter(p, 'LoadSessionData', 0);
parse(p, varargin{:});
sid = p.Results.sid;
annotFileName = p.Results.AnnotFile;
imgMetadataFileName = p.Results.imgMetadataFile;
refImgFileName = p.Results.refImgFile;
loadSessionData = p.Results.LoadSessionData;
wholeSession = [];

% Get sid if not provided
sid = str2double(regexp(sessionDataFile, '(?<=sid_).(?=_)', 'match'));

% Load imaging data session file
disp(['Loading ' sessionDataFile, '...'])
if loadSessionData
    load(fullfile(parentDir, sessionDataFile));
    sz = size(wholeSession);
    nPlanes = sz(3);
    nVolumes = sz(4);
    nTrials = sz(5);
else
    m = matfile(fullfile(parentDir, sessionDataFile));
    [~, ~, nPlanes, nVolumes, nTrials] = size(m, 'wholeSession');
end
disp('Session data loaded')

% Load imaging metadata file
imgMetadata = load(fullfile(parentDir, imgMetadataFileName)); % fields 'scanimageInfo', 'expDate' (scanimageInfo not present in older exps)
imgMetadata.sid = sid;

disp('Imaging metadata loaded')

% Add info from stim computer metadata files
stimDataFiles = dir(fullfile(parentDir, ['metadata*sid_', num2str(sid), '*.mat']));
imgMetadata.stimOnsetTimes = []; imgMetadata.stimDurs = []; imgMetadata.trialType = []; imgMetadata.outputData = [];
for iFile = 1:numel(stimDataFiles)
    
    % Load file    
    load(fullfile(parentDir, stimDataFiles(iFile).name)); % variable "metaData" with fields 'trialDuration', 'nTrials', 'stimTypes', 'sid', 'taskFile', 'outputData'
    
    % Add block data
    imgMetadata.blockData(iFile).nTrials = metaData.nTrials;
    imgMetadata.blockData(iFile).stimTypes = metaData.stimTypes;
    imgMetadata.blockData(iFile).taskFile = metaData.taskFile;
    imgMetadata.blockData(iFile).outputData = metaData.outputData;
    
    % Add info to overall session data
    imgMetadata.trialDuration = metaData.trialDuration;
    
    disp('Adding output data')
    % Add downsampled output data broken down into trials
    currOutput = metaData.outputData';
    sampPerTrial = size(currOutput, 2) / metaData.nTrials;
    rsOutput = reshape(currOutput, size(currOutput, 1), sampPerTrial, metaData.nTrials);   % --> [channel, sample, trial]
    disp(size(rsOutput))
    imgMetadata.outputData = cat(3, imgMetadata.outputData, permute(rsOutput(:,1:100:end,:), [2 1 3]));    % --> [sample, channel, trial]
    
    % Add trial type info
    disp('Adding trialType info');
    for iTrial = 1:metaData.nTrials
        currStim = metaData.stimTypes{iTrial};
        if strcmp(currStim, 'ClosedLoop')
            imgMetadata.stimOnsetTimes(end + 1) = nan;
            imgMetadata.stimDurs(end + 1) = nan;
            imgMetadata.trialType{end + 1} = 'ClosedLoop';
        else
            imgMetadata.stimOnsetTimes(end + 1) = str2double(regexp(currStim, '(?<=Onset_).*(?=_Dur)', 'match'));
            imgMetadata.stimDurs(end + 1) = str2double(regexp(currStim, '(?<=Dur_).*', 'match'));
            imgMetadata.trialType{end + 1} = regexp(currStim, '.*(?=_Onset)', 'match', 'once');
            if strcmp(imgMetadata.trialType{end}, 'NoOdor')
                imgMetadata.trialType{end} = 'NoStim'; % For backwards compatibility 
            end
        end
    end%iTrial
    
end%iFile

disp('Imaging data read')

% Load annotation data file
if ~isempty(annotFileName) && exist(fullfile(parentDir, annotFileName), 'file')
    disp(['Loading ' annotFileName, '...'])
    annotData = load(fullfile(parentDir, annotFileName)); % variables 'trialAnnotations', 'annotParams', 'ftData', 'flowArr', 'goodTrials', 'behaviorLabels', 'frameInfo'
    annotData.nFrames = annotData.frameInfo.nFrames;
    annotData.frameTimes = annotData.frameInfo.frameTimes;
    annotData.FRAME_RATE = annotData.frameInfo.FRAME_RATE; % This is the frame rate of the behavior video, not the GCaMP imaging
    disp([annotFileName, ' loaded'])
else
    disp('No behavioral annotation data loaded')
    annotData.trialAnnotations = [];
    annotData.annotParams = [];
    annotData.ftData = [];
    annotData.flowArr = [];
    annotData.goodTrials = [];
    annotData.behaviorLabels = [];
    annotData.frameInfo = [];
    annotData.nFrames = [];
    annotData.frameTimes = [];
    annotData.FRAME_RATE = [];
end

disp('annotation data loaded')

% Combine imaging data and annotation data into one structure
outputMetadata = setstructfields(imgMetadata, annotData);
outputMetadata.nTrials = nTrials;
outputMetadata.nPlanes = nPlanes;
outputMetadata.nVolumes = nVolumes;
outputMetadata.stimTypes = sort(unique(outputMetadata.trialType));

disp('Behavior data loaded')

% Extract metadata from scanimage info
if isfield(outputMetadata, 'scanimageInfo')
    outputMetadata.volumeRate = frameStringKeyLookup(outputMetadata.scanimageInfo.Software, 'SI.hRoiManager.scanVolumeRate');
    outputMetadata.laserPower = frameStringKeyLookup(outputMetadata.scanimageInfo.Software, 'SI.hBeams.powers');
    outputMetadata.scanFrameRate = frameStringKeyLookup(outputMetadata.scanimageInfo.Software, 'SI.hRoiManager.scanFrameRate');
else
    % For backwards compatibility
    outputMetadata.volumeRate = 6.44; % This is true for all older experiments
    outputMetadata.laserPower = [];
    outputMetadata.scanFrameRate = [];
end

% Separate trials by stim type
outputMetadata.stimSepTrials = [];
for iStim = 1:length(outputMetadata.stimTypes)
    outputMetadata.stimSepTrials.(outputMetadata.stimTypes{iStim}) = logical(cellfun(@(x) ...
        strcmp(x, outputMetadata.stimTypes{iStim}), outputMetadata.trialType));
end

% Match frame times to volumes if annotation data was provided (contains the video frame that
% most closely matches the time of the volume for each volume in the trial)
if ~isempty(outputMetadata.trialAnnotations)
    volTimes = (1:outputMetadata.nVolumes)' ./ outputMetadata.volumeRate;
    frameTimes = outputMetadata.frameTimes;
    for iVol = 1:outputMetadata.nVolumes
        [~, volFrames(iVol)] = min(abs(frameTimes - volTimes(iVol)));
        volFrames = volFrames';
    end
    outputMetadata.volFrames = volFrames;
else
    outputMetadata.volFrames = 1:outputMetadata.nVolumes;
end

disp('Loading reference images');

% Load a reference images file
if ~isempty(refImgFileName)
    refImages = load(fullfile(parentDir, refImgFileName));
    outputMetadata.refImg = refImages.refImages;
else
    disp('No reference images file selected')
    outputMetadata.refImg = [];
end

% Order fields alphabetically
outputMetadata = orderfields(outputMetadata);

catch ME
    write_to_log(ME.message, mfilename);
end

end