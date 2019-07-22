
%% LOAD IMAGING DATA AND BEHAVIORAL ANNOTATION DATA

expDir = uigetdir('D:\Dropbox (HMS)\2P Data\Imaging Data\', 'Select an imaging data directory');
savePath = fullfile(expDir, 'Analysis');

try
if expDir == 0
    % Skip loading if user clicked "Cancel"
    disp('Initialization cancelled')
else
    % Create analyis directory
    if ~isdir(savePath)
        mkdir(savePath)
    end
    
    % Load analysis metadata
    disp('Loading analysis metadata...')
    load(fullfile(expDir, 'analysisMetadata.mat')) % --> 'analysisMetadata' info struct
    infoStruct = analysisMetadata;

    % Load reference images file
    disp('Loading reference images...')
    if exist(fullfile(expDir, 'refImages_Reg.mat'), 'file')
        load(fullfile(expDir, 'refImages_Reg.mat')) % --> 'refImages', 'channelNum'
    else
        [refImgFile, refImgPath, ~] = uigetfile('*.mat', 'Select a reference images file', expDir);
        load(fullfile(regImgPath, refImgFile));
    end
    
    % Load PCA data
    disp('Loading PCA data...')
    dataFile = ['rigid_sid_', num2str(infoStruct.sid), '_sessionFile.mat'];
    if exist(fullfile(expDir, ['PCA_data_', dataFile ]), 'file')
        load(fullfile(expDir, ['PCA_data_', dataFile ])) % --> 'explained', 'pcaData' ([pc, plane], [y, x, pc, plane]
    else
        [pcaFile, pcaFilePath, ~] = uigetfile('*.mat', 'Select a PCA data file', expDir);
        if pcaFile
            load(fullfile(pcaFilePath, pcaFile)); % --> 'explained', 'pcaData' ([pc, plane], [y, x, pc, plane]
        end
    end
    
    % Load annotation type data
    disp('Loading annotation type data...')
    load(fullfile(expDir, 'annotationTypes.mat'))% --> 'annotationTypes', 'annotationTypeSummary'
    
    % Load FicTrac data
    if ~isfield(infoStruct, 'ftData')
        disp('Loading FicTrac data...')
        ftDir = fullfile('D:\Dropbox (HMS)\2P Data\Behavior Vids\', infoStruct.expDate, ...
            '_Movies\FicTracData');
        if isdir(ftDir)
            ftData = load_fictrac_data(infoStruct, 'Sid', infoStruct.sid, 'ParentDir', ftDir);
        else
            ftData = load_fictrac_data(infoStruct, 'Sid', infoStruct.sid);
        end
        infoStruct.ftData = ftData;
        analysisMetadata = infoStruct;
        save(fullfile(expDir, 'analysisMetadata.mat'), 'analysisMetadata', '-v7.3');
    else
        disp('FicTrac data already exists in infoStruct...skipping FicTrac loading process')
        ftData = infoStruct.ftData;
    end
    infoStruct.goodTrials(logical(ftData.resets)) = 0; % Omit any trials in which FicTrac reset
    
    % Load volume-averaged raw session data
    volAvgDataFile = fullfile(expDir, ['sid_', num2str(infoStruct.sid), ...
        '_volAvgSessionData.mat']);
    if exist(volAvgDataFile, 'file')
        load(volAvgDataFile) % "volAvgSessionData"
    end
    
    % Reload analysis metadata in case .goodTrials has been updated
    load(fullfile(expDir, 'analysisMetadata.mat')) % --> 'analysisMetadata' info struct
    infoStruct = analysisMetadata;
    infoStruct.refImg = refImages;
    disp('All data loaded')
    
    % ------- Load saved parameters for this experiment if they exist --------------------------
    p = load_plotting_params(expDir);
    if numel(fieldnames((p)))
        stimNames = p.stimNames; stimTrialGroups = p.stimTrialGroups;
        stimGroupNames = p.stimGroupNames; stimEpochs = p.stimEpochs;
        stimShadingColors = p.stimShadingColors; stimEpochNames = p.stimEpochNames;
        groupBounds = p.groupBounds; blockNames = p.blockNames; blockShading = p.blockShading;
    end
    
    % ------- Copy variables for convenience -------
    expDate = infoStruct.expDate;
    sid = infoStruct.sid;
    nPlanes = infoStruct.nPlanes;
    nVolumes = infoStruct.nVolumes;
    refImg = infoStruct.refImg;
    if ~isempty(infoStruct.nFrames)
        nFrames = infoStruct.nFrames;
    else
        nFrames = nVolumes;
    end
    nTrials = infoStruct.nTrials;
    nGoodTrials = sum(infoStruct.goodTrials);
    stimTypes = infoStruct.stimTypes;
    stimOnsetTimes = infoStruct.stimOnsetTimes;
    stimDurs = infoStruct.stimDurs;
    trialDuration = infoStruct.trialDuration;
    volumeRate = infoStruct.volumeRate;
    volFrames = infoStruct.volFrames;
    goodTrials = infoStruct.goodTrials;
    stimSepTrials = infoStruct.stimSepTrials; s = stimSepTrials;
    behaviorAnnotArr = annotationTypes{contains(annotationTypeSummary.AnnotationType, ...
        'move')}.frameAnnotArr;
    
    % Create hardcoded parameters
    FRAME_RATE = 25; % This is the frame rate for the behavior video
    MAX_INTENSITY = infoStruct.MAX_INTENSITY;
    if isempty(nFrames)
        nFrames = sum(trialDuration) * FRAME_RATE;
    end
    volTimes = (1:nVolumes)' ./ volumeRate;
    frameTimes = (1:nFrames)' ./ FRAME_RATE;
    
end%if
catch foldME; rethrow(foldME); end

%% LOAD BEHAVIOR DATA AND METADATA ONLY FOR OPTO STIM EXPERIMENT

clear variables; % to make sure variables don't get mixed up with an imaging experiment

sid = 0;
expDir = uigetdir('D:\Dropbox (HMS)\2P Data\Behavior Vids\', 'Select an experiment directory');
expDate = regexp(expDir, '201.*', 'match'); expDate = expDate{:};
savePath = fullfile(expDir, ['sid_', num2str(sid)]);

annotFileName = '_Movies\autoAnnotations.mat';

try
 
% Create analyis directory
if ~isdir(savePath)
    mkdir(savePath)
end
    
% ----------------  Load stim metadata -------------------------------------------------------------
infoStruct = [];
stimDataFiles = dir(fullfile(expDir, ['metadata*sid_', num2str(sid), '*.mat']));
infoStruct.stimOnsetTimes = []; infoStruct.stimDurs = []; infoStruct.trialType = []; infoStruct.outputData = [];
disp('Loading experiment metadata...')
for iFile = 1:numel(stimDataFiles)
    
    % Load file  
    load(fullfile(expDir, stimDataFiles(iFile).name)); % variable "metaData" with fields 'trialDuration', 'nTrials', 'stimTypes', 'sid', 'taskFile', 'outputData'
    
    % Add block data
    infoStruct.blockData(iFile).nTrials = metaData.nTrials;
    infoStruct.blockData(iFile).stimTypes = metaData.stimTypes;
    infoStruct.blockData(iFile).taskFile = metaData.taskFile;
    infoStruct.blockData(iFile).outputData = metaData.outputData;
    
    % Add info to overall session data
    infoStruct.trialDuration = metaData.trialDuration;
    infoStruct.expDate = expDate;
    
    % Add downsampled output data broken down into trials
    currOutput = metaData.outputData';
    sampPerTrial = size(currOutput, 2) / metaData.nTrials;
    rsOutput = reshape(currOutput, size(currOutput, 1), sampPerTrial, metaData.nTrials);   % --> [channel, sample, trial]
    disp(size(rsOutput))
    infoStruct.outputData = cat(3, infoStruct.outputData, permute(rsOutput(:,1:100:end,:), [2 1 3]));    % --> [sample, channel, trial]
    
    % Add trial type info
    for iTrial = 1:metaData.nTrials
        currStim = metaData.stimTypes{iTrial};
        infoStruct.stimOnsetTimes(end + 1) = str2double(regexp(currStim, '(?<=Onset_).*(?=_Dur)', 'match'));
        infoStruct.stimDurs(end + 1) = str2double(regexp(currStim, '(?<=Dur_).*', 'match'));
        infoStruct.trialType{end + 1} = regexp(currStim, '.*(?=_Onset)', 'match', 'once');
        if strcmp(infoStruct.trialType{end}, 'NoOdor')
            infoStruct.trialType{end} = 'NoStim'; % For backwards compatibility
        end
    end%iTrial
end%iFile
infoStruct.nTrials = size(infoStruct.outputData, 3);
infoStruct.stimTypes = sort(unique(infoStruct.trialType));
infoStruct.stimSepTrials = [];
for iStim = 1:length(infoStruct.stimTypes)
    infoStruct.stimSepTrials.(infoStruct.stimTypes{iStim}) = logical(cellfun(@(x) ...
        strcmp(x, infoStruct.stimTypes{iStim}), infoStruct.trialType));
end

% ----------------  Load autoAnnotation data -------------------------------------------------------
disp('Loading annotation data...')
annotData = load(fullfile(expDir, annotFileName)); % variables 'trialAnnotations', 'annotParams', 'ftData', 'flowArr', 'goodTrials', 'behaviorLabels', 'frameInfo'
annotData.nFrames = annotData.frameInfo.nFrames;
annotData.frameTimes = annotData.frameInfo.frameTimes;
annotData.FRAME_RATE = annotData.frameInfo.FRAME_RATE;
infoStruct = setstructfields(infoStruct, annotData);


% ----------------  Load FicTrac data --------------------------------------------------------------
disp('Loading FicTrac data...')
ftData = load_fictrac_data(infoStruct, 'sid', sid, 'ParentDir', fullfile(expDir, '\_Movies\FicTracData'));
infoStruct.ftData = ftData;
infoStruct.goodTrials(logical(ftData.resets)) = 0;

% ----------------  Load saved parameters for this experiment if they exist ------------------------
p = load_plotting_params(expDir);
if numel(fieldnames((p)))
    stimNames = p.stimNames; stimTrialGroups = p.stimTrialGroups; stimGroupNames = p.stimGroupNames;
    stimEpochs = p.stimEpochs; stimShadingColors = p.stimShadingColors; 
    stimEpochNames = p.stimEpochNames; groupBounds = p.groupBounds; blockNames = p.blockNames; 
    blockShading = p.blockShading;   
end

% ---------------- Create workspace vars -----------------------------------------------------------
infoStruct = orderfields(infoStruct);
nTrials = infoStruct.nTrials;
nFrames = infoStruct.nFrames;
stimTypes = infoStruct.stimTypes;
stimOnsetTimes = infoStruct.stimOnsetTimes;
stimDurs = infoStruct.stimDurs;
trialDuration = infoStruct.trialDuration;
goodTrials = infoStruct.goodTrials;
stimSepTrials = infoStruct.stimSepTrials; s = stimSepTrials;
behaviorAnnotArr = infoStruct.trialAnnotations;
FRAME_RATE = infoStruct.FRAME_RATE;
if ~isdir(savePath)
    mkdir(savePath)
end
catch foldME; rethrow(foldME); end

%% Plot raw fluorescence minimum for each trial
try
rawFlAvg = squeeze(min(min(volAvgSessionData, [], 1), [], 2)); % --> [plane, trial]
% rawFlAvg = squeeze(mean(mean(volAvgSessionData, 1), 2)); % --> [plane, trial]


flThresh = []; 
omitTrials = [];
% omitTrials = [1:10 15 17 22 25:29 32:35 39 80 82 91 101];
% omitTrials = find(rawFlAvg(5,:) > 50 | rawFlAvg(5,:) < 35)

if ~isempty(omitTrials)
    rawFlAvg(:,omitTrials) = [];
end
figure(1);clf;hold on
plot(rawFlAvg');
if ~isempty(flThresh)
   ylim([0 flThresh]) 
end

catch foldME; rethrow(foldME); end

%% SET UP PLOTTING PARAMETERS
try
stimNames = {'EtOH\_neat', 'OptoStim'}; % {'EtOH\_e-2'};%
stimTrialGroups = [s.OdorA + 2*s.OptoStim];%[ones(1, 30), 2 * ones(1, 60)];%[s.OdorA + 2 * s.OdorB]; % 
stimGroupNames = {'OdorA', 'OptoStim'};%{'OdorA', 'OdorB'}; %

stimEpochs = [10 15; 10 10.1];%[reshape(2:2:20, 2, 5)']% reshape(1:20, 2, 10)'% [11 12]; ;%[10 13];%[6 7;10 13];%[6 7; 8 11; 10 11];% [7 10 ; 10 13];%
stimShadingColors = {'red', 'green'};%repmat({'red'}, 1, size(stimEpochs, 1));%{'red'}; %{'red', 'green', 'red'}; % 
% stimShadingColors{3} = 'green';
stimEpochNames = {'Odor', 'OptoStim'};% {'Laser', 'Odor', 'Laser'}; %{'Laser', 'Odor'};%

%     groupBounds = [1:40:nTrials]; groupBounds(2:end) = groupBounds(2:end) - 1;
% groupBounds = [1, 40:60:nTrials-1]; groupBounds(2:end) = groupBounds(2:end) -1;
groupBounds = [0:10:nTrials-1];
groupBounds(1) = 0;

blockNames = {'Odor1', 'OptoStim1', 'Odor2', 'OptoStim2', 'Odor3', 'OptoStim3', 'Odor4', 'OptoStim4', 'Odor5', 'OptoStim5',}; %{'OdorOnly'};%{'Baseline', 'Photostim', 'OdorOnly'};%{'OdorOnly', 'BW', 'OdorOnly2', 'FW', 'OdorOnly3', 'BW2'}; %{'Odor Only', 'Photostim'}; %{'OdorOnly', 'FW', 'OdorOnly2', 'BW', 'OdorOnly3', 'FW2'};
blockShading = repmat({1,2}, 1, 5);%{2, [1 2], 2};%{2, [1 2], 2, [2 3], 2, [1 2]}; %{2, [2 3], 2, [1 2], 2, [2 3]};


stimGroupShading = {1, 2};

% Save plotting parameters for this experiment
overwrite = confirm_save(fullfile(expDir, 'plotting_params.mat'), 'DialogText', ...
    'This will overwrite an existing plotting parameter file...overwrite?');
if overwrite
    save(fullfile(expDir, 'plotting_params.mat'), 'stimNames', 'stimTrialGroups', 'stimGroupNames', ...
        'stimEpochs', 'stimShadingColors', 'stimEpochNames', 'groupBounds', 'blockNames', ...
        'blockShading');
    disp('Plotting parameters saved')
end

catch foldME; rethrow(foldME); end

%% =================================================================================================
%   BEHAVIOR SUMMARIES                                   
%%%=================================================================================================
     
    %% PLOT 2-D SUMMARY OF BEHAVIOR DATA ANNOTATIONS
try
saveDir = uigetdir(savePath, 'Select a save directory');
sepBlockStims = 0; sepGroupStims = 0;
% 
trialGroups = [goodTrials];
plotTitleSuffix = '';
fileNameSuffix = '_AllTrials';

% % % % % % % % % % % % % % % % plotTitleSuffix = ['  —  sid\_', num2str(sid)];

% % GROUP BY STIM TYPE
% trialGroups = stimTrialGroups; 
% trialGroups = trialGroups .* goodTrials;
% plotTitleSuffix = make_plotTitleSuffix(stimNames); %
% fileNameSuffix = make_fileNameSuffix(stimGroupNames);
% sepGroupStims = 1;

% % GROUP CHRONOLOGICALLY BY BLOCKS
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(groupBounds(end)+1:end) = iBound + 1;
% trialGroups = trialGroups .* goodTrials;
% plotTitleSuffix = '';
% fileNameSuffix = '_Blocks_Separated';
% if ~isempty(blockShading)
%     sepBlockStims = 1;
% end

% % GROUP BY ALTERNATING TRIALS
% trialGroups = zeros(1, nTrials);
% trialGroups(1:groupBounds(2)) = 1; % Baseline period
% trialGroups((groupBounds(2) + 1):2:end) = 2;
% trialGroups((groupBounds(2) + 2):2:end) = 3;
% trialGroups = trialGroups .* goodTrials;
% plotTitleSuffix = make_plotTitleSuffix({'Baseline', 'OdorOnly', 'Photostim'});
% fileNameSuffix = '_Alternating_Trials';
% if ~isempty(blockShading)
%     sepBlockStims = 1;  
% end


try

% Create plot titles
plotName = 'Behavior Annotation';
titleString = [regexprep(expDate, '_', '\\_'), '    ', [plotName, ' summary ', plotTitleSuffix]]; % regex to add escape characters
annotArr = behaviorAnnotArr;

% Create custom colormap
cMap = [rgb('Indigo'); ...
        rgb('Orange'); ...
        rgb('Green');
        rgb('Cyan'); ...
        ];
annotArr(end, end) = 4; % to keep the color mapping consistent

% Create figure
f = figure(7);clf
f.Position = [-1050 45 1020 950];
f.Color = [1 1 1];
ax = gca;
[~, ax, f] = plot_2D_summary(annotArr, FRAME_RATE, ...
    'plotAxes', ax, ...
    'trialGroups', trialGroups, ...
    'titleStr', titleString, ...
    'colormap', cMap ...
    );

ax.FontSize = 14;
ax.Title.FontSize = 12;
ax.XLabel.FontSize = 14;

% Plot stim times
hold on
if sepBlockStims
    
    % Calculate y-axis ranges for each block
    blockStarts = 1;
    blockEnds = [];
    for iBlock = 1:numel(unique(trialGroups(trialGroups > 0)))-1
        blockTrials = sum(trialGroups == iBlock);
        blockStarts(iBlock + 1) = blockStarts(iBlock) + 4 + blockTrials;
        blockEnds(iBlock) = blockStarts(iBlock) + blockTrials - 1;
    end
    blockEnds(end + 1) = numel(trialGroups) + (4 * numel(unique(trialGroups(trialGroups > 0))));
    blockRanges = [blockStarts - 1; blockEnds + 1]';
    
    % Plot stim start and end for each block
    for iBlock = 1:numel(unique(trialGroups(trialGroups > 0)))
        currBlockShading = stimEpochs(blockShading{iBlock}, :);
        nStimEpochs = size(currBlockShading, 1);
        currShadeColors = stimShadingColors(blockShading{iBlock});
        draw_stim_lines(currBlockShading, currShadeColors, ...
                        'plotAxes', ax, 'yLims', blockRanges(iBlock, :));      
    end
    
elseif sepGroupStims
    
    % Calculate y-axis ranges for each stim type
    blockStarts = 1;
    blockEnds = [];
    for iBlock = 1:numel(unique(trialGroups(trialGroups > 0)))-1
        blockTrials = sum(trialGroups == iBlock);
        blockStarts(iBlock + 1) = blockStarts(iBlock) + 4 + blockTrials;
        blockEnds(iBlock) = blockStarts(iBlock) + blockTrials - 1;
    end
    blockEnds(end + 1) = numel(trialGroups) + (4 * numel(unique(trialGroups(trialGroups > 0))));
    blockRanges = [blockStarts - 1; blockEnds + 1]';
    
    % Plot stim start and end for each block
    for iBlock = 1:numel(unique(trialGroups(trialGroups > 0)))
        currBlockShading = stimEpochs(stimGroupShading{iBlock}, :);
        nStimEpochs = size(currBlockShading, 1);
        currShadeColors = stimShadingColors(stimGroupShading{iBlock});
        draw_stim_lines(currBlockShading, currShadeColors, ...
                        'plotAxes', ax, 'yLims', blockRanges(iBlock, :));      
    end
    
else
    % Draw one set of lines for the entire plot
    draw_stim_lines(stimEpochs, stimShadingColors, 'plotAxes', ax);
end

% Save figure
if saveDir
    fileName = ['2D_Behavior_Annotation_Summary', regexprep(expDate, {'_', 'exp'}, {'', '_'}), ...
                fileNameSuffix];
    save_figure(f, saveDir, fileName);
end%if

catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end
    %% PLOT 1-D SUMMARY OF BEHAVIOR DATA ANNOTATIONS
try
%----- Plot 1D trial-averaged movement data -----
s = stimSepTrials; 
actionNames = {'NA', 'IsoMovement', 'Locomotion'};
saveDir = uigetdir(savePath, 'Select a save directory');
actionNum = [3]; % locomotionLabel = 3; noActionLabel = 0; isoMovementLabel = 1;
actionName = actionNames{actionNum};
figTitle = regexprep([expDate, '  —  Fly ', actionName, ' throughout trial'], '_', '\\_');
cm = [];

% % ALL TRIALS
% trialGroups = [goodTrials];
% fileNameSuffix = ['_AllTrials_', actionName];
% groupNames = {'All trials'};
% 
% figTitle = [figTitle, ['  —  sid\_', num2str(sid)]];

% GROUP BY STIM TYPE
trialGroups = stimTrialGroups .* goodTrials;
fileNameSuffix = [make_fileNameSuffix(stimGroupNames), '_', actionName]; 
groupNames = stimNames;
% % 
% % GROUP BY EARLY/LATE
% groupNames = [];
% groupBounds = [1, 30, 60];
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
%    groupNames{iBound} = ['Trials ', num2str(groupBounds(iBound)), '-', num2str(groupBounds(iBound + 1))];
% end
% trialGroups(groupBounds(end):end) = numel(groupBounds);
% groupNames{end + 1} = ['Trials ', num2str(groupBounds(end)), '-', num2str(nTrials)];
% fileNameSuffix = '_EarlyVsLateTrials';
% 
% % PLOT AND COLOR EACH BLOCK SEPARATELY
% groupNames = blockNames;
% trialGroups = zeros(1, nTrials);
% trialGroups = trialGroups .* goodTrials;
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(groupBounds(end)+1:end) = numel(groupBounds);
% fileNameSuffix = '_Blocks_Separated';

% % GROUP BY BLOCK TYPE 
% groupNames = blockNames(1:2); 
% trialGroups = zeros(1, nTrials);
% trialGroups = trialGroups .* goodTrials;
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(logical(mod(trialGroups, 2))) = 1;
% trialGroups(~logical(mod(trialGroups, 2))) = 2;
% fileNameSuffix = '_Block_Type_Groups';

% % 
% % GROUP BY ALTERNATING TRIALS (WITH BASELINE)
% trialGroups = zeros(1, nTrials);
% trialGroups(1:groupBounds(2)) = 1; % Baseline period
% trialGroups((groupBounds(2) + 1):2:end) = 2;
% trialGroups((groupBounds(2) + 2):2:end) = 3;
% trialGroups = trialGroups .* goodTrials;
% fileNameSuffix = '_Alternating_Trials';
% groupNames = blockNames;
% if ~isempty(blockShading)
%     sepBlockStims = 1;  
% end
% 
%     % Drop the baseline block
%     groupNames = blockNames(2:3);
%     trialGroups(1:40) = 0;
%     trialGroups = trialGroups - 1;
%     trialGroups(trialGroups < 0) = 0;
%     cm = [rgb('red'); rgb('green'); rgb('blue')];
%     fileNameSuffix = '_Alternating_Trials_No_Baseline';

try

% Create array of annotation data
f = figure(2); clf; hold on
f.Position = [-1600 300 1600 500];
f.Color = [1 1 1];

if isempty(trialGroups)
    
    % Plot summed movement data
    annotArrSum = sum(ismember(behaviorAnnotArr, actionNum), 1) ./ nTrials;
    ax = gca();
    ax.FontSize = 14;
    plot_behavior_summary_1D(infoStruct, annotArrSum(2:end-1), ax, figTitle);
    
else
    annotArrSum = [];
    yLimsAll = [];
    ax = [];
    nGroups = length(unique(trialGroups(trialGroups ~= 0)));
    if isempty(cm)
        cm = parula(nGroups);
        cm = [rgb('blue'); rgb('red'); rgb('green'); rgb('magenta'); rgb('cyan'); rgb('gold'); ...
                rgb('lime'); rgb('black'); rgb('maroon'); rgb('grey')];
    end
    for iGroup = 1:nGroups
        ax = gca();
        ax.FontSize = 14;
        colormap(jet(nGroups))
        annotArrSum = sum(ismember(behaviorAnnotArr(trialGroups == iGroup, :), actionNum), 1) ./ sum(trialGroups == iGroup);
        [plt, ~, ~] = plot_behavior_summary_1D(infoStruct, annotArrSum, 'PlotAxes', ax, 'LineColor', cm(iGroup, :));
        plt.LineWidth = 2;
        if iGroup ~= length(unique(trialGroups))
            xlabel('');
        end
    end%iGroup
    
    ax.XLim = [10 nFrames-10]; % to improve plot appearance
    ax.YLim = [0 1];
    suptitle(figTitle);
end

% Add shading during stimulus presentations
yL = ylim();
for iStim = 1:size(stimEpochs, 1)
    stimOnset = stimEpochs(iStim, 1) * FRAME_RATE;
    stimOffset = stimEpochs(iStim, 2) * FRAME_RATE;
    plot_stim_shading([stimOnset, stimOffset], 'Color', rgb(stimShadingColors{iStim}))
end
legend(cat(2, groupNames, unique(stimEpochNames)), 'FontSize', 14, 'Location', 'SE', 'AutoUpdate', 'off')

% Save figure
if saveDir
    fileName = regexprep(['1D_Behavior_Annotation_Summary_', ...
            regexprep(expDate, {'_', 'exp'}, {'', '-'}), fileNameSuffix], '_', '\_');
    save_figure(f, saveDir, fileName);
end%if

catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end
    %% PLOT 2-D SUMMARY OF FICTRAC DATA
try
saveDir = uigetdir(savePath, 'Select a save directory');
fontSize = 14;
sepBlockStims = 0; sepGroupStims = 0;

% ftVarName = 'moveSpeed'; % 'moveSpeed', 'fwSpeed', 'yawSpeed', 'yawVel'
ftVarName =  'fwSpeed'%;%'yawVel' %              'yawSpeed'%  'fwSpeed'% 
sdCap = 300;
smWin = 1;
cmName = @parula;
figTitle = [regexprep(expDate, '_', '\\_'), '  —  FicTrac ', ftVarName];

% % % % % 
% ALL TRIALS
trialGroups = [];
fileNameSuffix = ['_AllTrials'];
plotTitleSuffix = '';
% 
% plotTitleSuffix = ['  —  sid\_', num2str(sid)];
% % 
% % GROUP CHRONOLOGICALLY BY BLOCKS
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(groupBounds(end)+1:end) = iBound + 1;
% trialGroups = trialGroups .* goodTrials;
% plotTitleSuffix = '';
% fileNameSuffix = '_Blocks_Separated';
% if ~isempty(blockShading)
%     sepBlockStims = 1;
% end
% 
% % GROUP BY STIM TYPE
% trialGroups = stimTrialGroups; 
% trialGroups = trialGroups .* goodTrials;
% plotTitleSuffix = make_plotTitleSuffix(stimNames); %
% fileNameSuffix = make_fileNameSuffix(stimGroupNames);
% sepGroupStims = 1;

% % % 
% % GROUP BY ALTERNATING TRIALS
% trialGroups = zeros(1, nTrials);
% trialGroups(1:groupBounds(2)) = 1; % Baseline period
% trialGroups((groupBounds(2) + 1):2:end) = 2;
% trialGroups((groupBounds(2) + 2):2:end) = 3;
% trialGroups = trialGroups .* goodTrials;
% plotTitleSuffix = make_plotTitleSuffix({'Baseline', 'Photostim', 'OdorOnly'});
% fileNameSuffix = '_Alternating_Trials';
% if ~isempty(blockShading)
%     sepBlockStims = 1;  
% end

try
    
% Extract FicTrac data
if strcmp(ftVarName, 'yawSpeed')
    rawData = infoStruct.ftData.(ftVarName);          % --> [frame, trial]
    rawData = rawData';                                     % --> [trial, frame]
    plotData = abs(rad2deg(rawData .* FRAME_RATE));    	% --> [trial, frame] (deg/sec)
    figTitle = [figTitle, ' (deg/sec)'];
elseif strcmp(ftVarName, 'yawVel')
    rawData = infoStruct.ftData.yawSpeed;          % --> [frame, trial]
    rawData = rawData';                                     % --> [trial, frame]
    plotData = (rad2deg(rawData .* FRAME_RATE));    	% --> [trial, frame] (deg/sec)
    figTitle = [figTitle, ' (deg/sec)'];
else
    rawData = infoStruct.ftData.(ftVarName);          % --> [frame, trial]
    rawData = rawData';
    plotData = rawData .* FRAME_RATE .* 4.5;            % --> [trial, frame] (mm/sec)
    figTitle = [figTitle, ' (mm/sec)'];
end

% Cap values at n SD above mean
capVal = mean(plotData(:), 'omitnan') + (sdCap * std(plotData(:), 'omitnan'));
h = figure(2); clf; h.Position = [-1900 100 800 800];
subplot(211); hist(as_vector(smoothdata(plotData, 2, 'gaussian', smWin)), 100);
plotData(plotData > capVal) = capVal;
subplot(212);  hist(as_vector(smoothdata(plotData, 2, 'gaussian', smWin)), 100)

% Smooth data
smPlotData = smoothdata(plotData, 2, 'gaussian', smWin);

% Create colormap
cm = cmName(numel(unique(smPlotData)));

% Plot data
titleStr = [figTitle, plotTitleSuffix];
fileName = ['2D_FicTrac_', ftVarName '_Summary', fileNameSuffix, '_', ...
            regexprep(expDate, {'_', 'exp'}, {'', '_'})];
[~, ax, f] = plot_2D_summary(smPlotData, FRAME_RATE, ...
                'trialGroups', trialGroups, ...
                'titleStr', titleStr, ...
                'colormap', cm ...
                );
if sepBlockStims
   ax.YTick = []; 
end            
ax.Title.FontSize = fontSize;
ax.FontSize = fontSize;
% Update colormap if necessary
if strcmp(ftVarName, 'yawVel')
    colormap(gca, bluewhitered)
end

colorbar
hold on

% Plot stim times
hold on
if sepBlockStims
    
    % Calculate y-axis ranges for each block
    blockStarts = 1;
    blockEnds = [];
    for iBlock = 1:numel(unique(trialGroups(trialGroups > 0)))-1
        blockTrials = sum(trialGroups == iBlock);
        blockStarts(iBlock + 1) = blockStarts(iBlock) + 4 + blockTrials;
        blockEnds(iBlock) = blockStarts(iBlock) + blockTrials - 1;
    end
    blockEnds(end + 1) = numel(trialGroups) + (4 * numel(unique(trialGroups(trialGroups > 0))));
    blockRanges = [blockStarts - 1; blockEnds + 1]';
    
    % Plot stim start and end for each block
    for iBlock = 1:numel(unique(trialGroups(trialGroups > 0)))
        currBlockShading = stimEpochs(blockShading{iBlock}, :);
        nStimEpochs = size(currBlockShading, 1);
        currShadeColors = stimShadingColors(blockShading{iBlock});
        draw_stim_lines(currBlockShading, currShadeColors, ...
            'plotAxes', ax, 'yLims', blockRanges(iBlock, :));
    end
elseif sepGroupStims
    
    % Calculate y-axis ranges for each stim type
    blockStarts = 1;
    blockEnds = [];
    for iBlock = 1:numel(unique(trialGroups(trialGroups > 0)))-1
        blockTrials = sum(trialGroups == iBlock);
        blockStarts(iBlock + 1) = blockStarts(iBlock) + 4 + blockTrials;
        blockEnds(iBlock) = blockStarts(iBlock) + blockTrials - 1;
    end
    blockEnds(end + 1) = numel(trialGroups) + (4 * numel(unique(trialGroups(trialGroups > 0))));
    blockRanges = [blockStarts - 1; blockEnds + 1]';
    
    % Plot stim start and end for each block
    for iBlock = 1:numel(unique(trialGroups(trialGroups > 0)))
        currBlockShading = stimEpochs(stimGroupShading{iBlock}, :);
        nStimEpochs = size(currBlockShading, 1);
        currShadeColors = stimShadingColors(stimGroupShading{iBlock});
        draw_stim_lines(currBlockShading, currShadeColors, ...
            'plotAxes', ax, 'yLims', blockRanges(iBlock, :));
    end
    
else
    % Draw one set of lines for the entire plot
    draw_stim_lines(stimEpochs, stimShadingColors, 'plotAxes', ax);
end

% Save figure
if saveDir
    fileName = ['2D_FicTrac_', ftVarName, '_Summary_', ...
                regexprep(expDate, {'_', 'exp'}, {'', '-'}), fileNameSuffix, '_', ];
    save_figure(f, saveDir, fileName);
end

catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end
    %% PLOT 1-D SUMMARY OF FICTRAC DATA
try
    
saveDir = uigetdir(savePath, 'Select a save directory');

includeQuiescence = 1;
if ~includeQuiescence
    fileNameSuffix = 'NoQuiescence_';
else
    fileNameSuffix = '';
end
figTitle = [expDate, '  —  Trial-Averaged FicTrac data'];
trialGroups = goodTrials';
smWin = 9;
cm = [];
% 
% 
% % ALL TRIALS
% trialGroups = [goodTrials];
% fileNameSuffix = [fileNameSuffix, 'AllTrials'];
% groupNames = {'All trials'};
% 
% figTitle = [figTitle, '  —  sid_', num2str(sid)];

% % % 
% % GROUP BY STIM TYPE
% trialGroups = stimTrialGroups .* goodTrials; 
% fileNameSuffix = [fileNameSuffix, '_StimTypeComparison']; 
% groupNames = stimNames;
% % 
% % % % 

% % GROUP BY EARLY/LATE
% binNames = [];
% binBounds = [1, 30, 70];
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(binBounds)-1
%     trialGroups(binBounds(iBound):binBounds(iBound + 1)) = iBound;
%     binNames{iBound} = ['Trials ', num2str(binBounds(iBound)), '-', num2str(binBounds(iBound + 1))];
% end
% trialGroups(binBounds(end):end) = numel(binBounds);
% trialGroups = trialGroups .* goodTrials;
% binNames{end + 1} = ['Trials ', num2str(binBounds(end)), '-', num2str(nTrials)];
% fileNameSuffix = [fileNameSuffix, 'EarlyVsLateTrials'];
% % 
% PLOT AND COLOR EACH BLOCK SEPARATELY
groupNames = blockNames;
trialGroups = zeros(1, nTrials);
for iBound = 1:numel(groupBounds)-1
   trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
end
trialGroups(groupBounds(end)+1:end) = numel(groupBounds);
trialGroups = trialGroups .* goodTrials;
fileNameSuffix = [fileNameSuffix, 'Blocks_Separated'];

% 
% trialGroups = trialGroups - 1;
% trialGroups(1:20) = 0;
% groupNames = groupNames(2:4);

% % 
% % GROUP BY BLOCK TYPE 
% groupNames = blockNames(1:2); 
% trialGroups = zeros(1, nTrials);
% trialGroups = trialGroups .* goodTrials;
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(logical(mod(trialGroups, 2))) = 1;
% trialGroups(~logical(mod(trialGroups, 2))) = 2;
% fileNameSuffix = [fileNameSuffix, 'Block_Type_Groups'];

% % GROUP BY ALTERNATING TRIALS (WITH BASELINE)
% trialGroups = zeros(1, nTrials);
% trialGroups(1:groupBounds(2)) = 1; % Baseline period
% trialGroups((groupBounds(2) + 1):2:end) = 2;
% trialGroups((groupBounds(2) + 2):2:end) = 3;
% trialGroups = trialGroups .* goodTrials;
% fileNameSuffix = [fileNameSuffix, 'Alternating_Trials'];
% groupNames = blockNames;
% if ~isempty(blockShading)
%     sepBlockStims = 1;  
% end
% % % % % 
%     % Drop the baseline block
%     groupNames = blockNames(2:3);
%     trialGroups(1:40) = 0;
%     trialGroups = trialGroups - 1;
%     trialGroups(trialGroups < 0) = 0;
%     cm = [rgb('red'); rgb('green'); rgb('blue')];
%     fileNameSuffix = [fileNameSuffix, '_No_Baseline'];
% %     
% trialGroups(100:end) = 0;

try
    
% Extract relevant data
xTickFR = [0:1/trialDuration:1] * (trialDuration * FRAME_RATE);
xTickLabels = [0:1/trialDuration:1] * trialDuration;
mmSpeedData = ftData.moveSpeed * FRAME_RATE * 4.5;  % --> [frame, trial] (mm/sec)
dHD = abs(rad2deg(ftData.yawSpeed * FRAME_RATE));        % --> [frame, trial] (deg/sec)
fwSpeed = (rad2deg(ftData.yawSpeed * FRAME_RATE));        % --> [frame, trial] (deg/sec)%ftData.fwSpeed * FRAME_RATE * 4.5;        % --> [frame, trial  (mm/sec)
yawVel = rad2deg(ftData.yawSpeed * FRAME_RATE);        % --> [frame, trial] (deg/sec)
nFrames = size(mmSpeedData, 1);

% Create figurenum
f = figure(3); clf; hold on
f.Position = [-1100 50 900 930];
f.Color = [1 1 1];

% Create axes
M = 0.02;
P = 0.00;
axMoveSpeed = subaxis(3,1,1, 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.05); hold on
axYawVel = subaxis(3,1,2, 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.06); hold on
axYawSpeed = subaxis(3,1,3, 'S', 0, 'M', M, 'PB', 0.06, 'PL', 0.06); hold on

% Plot data
nGroups = length(unique(trialGroups(trialGroups ~= 0)));
if isempty(cm)
    cm = parula(nGroups);
    cm = [rgb('blue'); rgb('red'); rgb('green'); rgb('magenta'); rgb('cyan'); rgb('gold'); ...
            rgb('lime');rgb('black');rgb('maroon');rgb('grey')];
end
for iGroup = 1:nGroups
    
    % Calculate mean values for current group
    currXYSpeed = mmSpeedData(:, trialGroups==iGroup);
    currYawSpeed = dHD(:, trialGroups == iGroup);
    currYawVel = yawVel(:, trialGroups == iGroup);
    
    if ~includeQuiescence
        currAnnotData = behaviorAnnotArr';
        currAnnotData(:, trialGroups ~= iGroup) = [];
        currXYSpeed(currAnnotData == 0) = nan;
        currYawSpeed(currAnnotData == 0) = nan;
        currYawVel(currAnnotData == 0) = nan;
    end
    
    % Omit outliers
    outlierCalc = @(x) mean(x) + 4 * std(x);
    currXYSpeed(currXYSpeed >= outlierCalc(mmSpeedData(:))) = nan;
    currYawSpeed(currYawSpeed >= outlierCalc(dHD(:))) = nan;
    currYawVel(currYawVel >= outlierCalc(yawVel(:))) = nan;

    meanSpeed = smoothdata(mean(currXYSpeed, 2, 'omitnan'), 'gaussian', smWin);
    meanYawSpeed = smoothdata(mean(currYawSpeed, 2, 'omitnan'), 'gaussian', smWin);
    meanYawVel = smoothdata(mean(currYawVel, 2, 'omitnan'), 'gaussian', smWin);

    % XY speed plot
    axes(axMoveSpeed)
    plot(meanSpeed, 'linewidth', 2, 'color', cm(iGroup, :));
    
    % Directional yaw velocity plot
    axes(axYawVel)
    plot(meanYawVel, 'linewidth', 2, 'color', cm(iGroup,:));
        
    % Yaw speed plot
    axes(axYawSpeed)
    plot(meanYawSpeed, 'linewidth', 2, 'color', cm(iGroup,:));
    axYawSpeed.XLabel.String = 'Time (sec)';

end%iGroup


% Format axes
allAxes = {axMoveSpeed, axYawVel, axYawSpeed};
yAxisLabels = {'XY Speed (mm/sec)', 'Yaw Vel (CW = +)', 'Yaw Speed (deg/sec)'};
for iAxis = 1:numel(allAxes)
    currAxis = allAxes{iAxis};
    currLabel = yAxisLabels{iAxis};
    
    currAxis.XTick = xTickFR;
    currAxis.XTickLabel = xTickLabels;
    currAxis.YLabel.String = currLabel;
    currAxis.FontSize = 14;
    legend(currAxis, groupNames, 'FontSize', 12, 'Location', 'NW', 'AutoUpdate', 'off')
    currAxis.XLim = [9 nFrames-5]; % to improve plot appearance
    if ~isempty(stimEpochs)
        for iStim = 1:size(stimEpochs, 1)
            stimOnset = stimEpochs(iStim, 1) * FRAME_RATE;
            stimOffset = stimEpochs(iStim, 2) * FRAME_RATE;
            plot_stim_shading([stimOnset, stimOffset], 'Color', rgb(stimShadingColors{iStim}), ...
                    'Axes', currAxis);
        end
    end
end

suptitle(regexprep([figTitle, '  —  ', fileNameSuffix], '_', '\\_'));

% Save figure
if saveDir
    fileName = regexprep(['1D_FicTrac_Summary_', regexprep(expDate, {'_', 'exp'}, {'', '1'}), ...
                '_', fileNameSuffix], '_', '\_');
    save_figure(f, saveDir, fileName);
end%if

catch foldME; rethrow(foldME); end

catch foldME; rethrow(foldME); end
 %% PLOT OVERLAID 2D MOVEMENT DATA
try
startTime = 15;
plotLen = 1;
limScalar = 0.8;
alpha = 0.4;
showMean = 1;
saveFig = 0;

s = stimSepTrials;
xTickFR = [0:1/trialDuration:1] * (trialDuration * FRAME_RATE);
xTickLabels = [0:1/trialDuration:1] * trialDuration;
surfPlot = 0;
mmXY = ftData.intXY * 4.5; % --> [frame, var, trial] (mm)
HD = ftData.intHD;    % --> [frame, trial]
nFrames = size(mmXY, 1);
startFrame = startTime * FRAME_RATE;
endFrame = (startTime + plotLen) * FRAME_RATE;
trialGroups = ones(size(mmXY, 3), 1) .* goodTrials';
cm = [rgb('blue'); rgb('red'); rgb('green'); rgb('magenta'); rgb('cyan'); rgb('gold'); rgb('lime');rgb('black')];

% % PLOT ALL TRIALS COLORED BY TIME
% surfPlot = 1;
% cm = jet(endFrame - startFrame + 1);
% fileNameSuffix = '_chronological';
% % 
% % GROUP BY STIM TYPE
% trialGroups =  [s.OdorA + 2 * s.OdorB] .* goodTrials;
% fileNameSuffix = '_OdorAvsOdorBvsNoStim'; 
% groupNames = stimNames;

% PLOT AND COLOR EACH BLOCK SEPARATELY
groupNames = blockNames;
trialGroups = zeros(1, nTrials);
for iBound = 1:numel(groupBounds)-1
   trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
end
trialGroups(groupBounds(end)+1:end) = numel(groupBounds);
trialGroups = trialGroups .* goodTrials;
fileNameSuffix = '_Blocks_Separated';

% % GROUP BY ALTERNATING TRIALS
% trialGroups = zeros(1, nTrials);
% trialGroups(1:groupBounds(2)) = 1; % Baseline period
% trialGroups((groupBounds(2) + 1):2:end) = 2;
% trialGroups((groupBounds(2) + 2):2:end) = 3;
% trialGroups = trialGroups .* goodTrials;
% fileNameSuffix = '_Alternating_Trials';
% groupNames = blockNames;
% if ~isempty(blockShading)
%     sepBlockStims = 1;  
% end
% 
%     % Drop the baseline block
%     groupNames = blockNames(2:3);
%     trialGroups(1:40) = 0;
%     trialGroups = trialGroups - 1;
%     trialGroups(trialGroups < 0) = 0;
%     cm = [rgb('red'); rgb('green'); rgb('blue')];

% trialGroups(100:end) = 0;

try 
%%% OVERLAY 2D MOVEMENT DATA

% Create and format figure
f = figure(13); clf; hold on; 
f.Color = [1 1 1];
f.Position = [-1035 50 1000 800];
ax = gca();
ax.FontSize = 20;
title(regexprep([expDate, ' -  Time Window: ', num2str(startTime), '-', num2str(startTime + plotLen), ' sec'], '_', '\\_'))

% Process movement data
maxPos = 0;
plotX = []; plotY = [];
mmXY = ftData.intXY * 4.5;
nGroups = numel(unique(trialGroups(trialGroups > 0)));
for iTrial = 1:size(mmXY, 3)%[16 80 28]
    
    currXY = mmXY(:, :, iTrial);
    smoothX = smoothdata(currXY(:, 1), 'gaussian', 7);
    smoothY = smoothdata(currXY(:, 2), 'gaussian', 7);
    
    % Rotate so fly is facing upwards at starting point
    smoothHD = mod(smoothdata(unwrap(HD(:, iTrial)), 'gaussian', 11), (2*pi));
    theta = -smoothHD(startFrame) + pi/2;
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    smoothXYRot = R * [smoothX'; smoothY'];
    xData = smoothXYRot(1, startFrame:endFrame);
    yData = smoothXYRot(2, startFrame:endFrame);
    
    % Center fly at (0,0)
    xData = xData - xData(1);
    yData = yData - yData(1);
    
    plotX(:, iTrial) = xData; % --> [frame, trial]
    plotY(:,iTrial) = yData;  % --> [frame, trial]

    % For setting axis limits later
    for iSec = 1:(trialDuration - plotLen - 1) 
        currStart = FRAME_RATE * iSec;
        currEnd = currStart + (FRAME_RATE * plotLen);
        currPosData = currXY(currStart:currEnd, :);
        currPosData(:,1) = currPosData(:,1) - currPosData(1,1);
        currPosData(:,2) = currPosData(:,2) - currPosData(1,2);
        maxPos = max([maxPos; currPosData(:)]);
    end    
end%iTrial

% Plot traces
if surfPlot
    for iTrial = 1:nTrials
        if goodTrials(iTrial) 
            % Plot a color gradient for each line
            currX = plotX(:,iTrial);
            currY = plotY(:,iTrial);
            cData = permute(repmat(cm, 1, 1, 2), [1 3 2]);
            surface('XData', [currX currX], ...
                'YData', [currY currY], ...
                'ZData', zeros(numel(currX), 2), ...
                'CData', cData, ...
                'FaceColor', 'none', 'EdgeColor', 'interp', 'marker', 'none');
        end%if
    end%iTrial
else
    % Plot colored by trial groups
    legendPlots = [0 0 0 0 0 0 0 0 0 0 0 0]; legendObj = [];
    for iGroup = 1:nGroups
        for iTrial = 1:nTrials
            if trialGroups(iTrial) == iGroup
                
                plt = plot(plotX(:,iTrial), plotY(:,iTrial), 'color', [cm(iGroup, :), alpha], 'linewidth', 1);
                
                % Save one plot line from each group to use in legend
                if ~legendPlots(trialGroups(iTrial))
                    legendObj(trialGroups(iTrial)) = plt;
                    legendPlots(trialGroups(iTrial)) = 1;
                end
            end%if
        end%iTrial
        
        if showMean
            for iGroup = 1:nGroups
                % Calculate mean ending vector in the current group
                endX = plotX(end, trialGroups == iGroup)';
                endY = plotY(end, trialGroups == iGroup)';
                [endThetas, endRhos] = cart2pol(endX, endY);
                mu = circ_mean(endThetas, endRhos);
                [meanX, meanY] = pol2cart(mu, mean(endRhos));
                plot([0 meanX], [0 meanY], 'color', cm(iGroup, :), 'linewidth', 3)
                altX = mean(endX);
                altY = mean(endY);
                plot([0 altX], [0 altY], 'color', cm(iGroup, :), 'linewidth', 3)
            end%iGroup
        end
    end%iGroup
    legend(legendObj, groupNames, 'FontSize', 14, 'Location', 'best');
end%if

% Set axis limits
lims = 1.05 * maxPos * limScalar;
axis equal
xlabel('2D movement (mm)');
ylabel('2D movement (mm)');
xlim([-lims lims])
ylim([-lims lims])
tightfig;
ax.FontSize = 14;


% Save figure if necessary
if saveFig
    saveDir = uigetdir(savePath, 'Select a save directory');
    if saveDir
        % Create filename
        fileName = regexprep(['Movement_Overlay_', regexprep(expDate, {'_', 'exp'}, {'', '_'}), ...
                '_', num2str(startTime), '-', num2str(startTime + plotLen)...
                '_sec', fileNameSuffix], '_', '\_');
        save_figure(f, saveDir, fileName);
    end
end
catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end

%% =================================================================================================
%   LOAD PROCESSED EVENT DATA 
%% =================================================================================================
try
% Load processed event data
[eventDataFile, eventDataPath, ~] = uigetfile('*.mat', 'Select an event data file', expDir);
if eventDataFile ~= 0
    load(fullfile(eventDataPath, eventDataFile))
    %===================================================================================================
    % --> 'alignEventSummary', 'filterEventSummary', 'primaryEventNames', 'eventLists',
    %     'nEventTypes', 'condNames', 'onsetFilterVecs', 'offsetFilterVecs', 'analysisWindows'
    %
    %
    %                                                 |------------event------------|
    %       <--------fW(1)--------><----aW(1)---->[alignVol]<----aW(2)----><----fW(2)---->
    %       [------------------fD(1)-----------------][-------fD(2)-------][----fD(3)----]
    %
    %===================================================================================================
    clear eventDataPath
    disp(primaryEventNames)
    disp(alignEventSummary)
    disp(filterEventSummary)
    
    % Load dF/F data for filtered alignment events
    load(fullfile(expDir, ['CombDffAvg_', eventDataFile])) % --> 'combinedDffAvg', 'allCondSummaries', 'allCondNames', 'combFilterVecs'
else
    disp('No event data file selected - loading cancelled')
end
catch foldME; rethrow(foldME); end
%% =================================================================================================
%   ODOR/SOUND STIM HEATMAPS                                
%%%=================================================================================================
try
    % Show summary again
    odorEventName = 'OdorA';
    eventInd = cellfun(@strcmp, primaryEventNames, repmat({odorEventName}, 1, numel(primaryEventNames)));
    currSummary = allCondSummaries{eventInd};
    currCondNames = allCondNames{eventInd};
    disp(currSummary)
    
    % Calculate absolute max dF/F value across all planes and stim types
    currConds = [1 4 2 5];
    makeVid = 1;
    sigma = [0.6];   
    rangeType = 'Max';
    rangeScalar = 0.4;
    saveDir = [];
    fileName = [odorEventName, '_Response_Heatmaps_', regexprep(expDate, {'_', 'exp'}, {'', '-'})];

    plotTitles = [];
    for iCond = currConds
        plotTitles{iCond} = [regexprep(currCondNames{iCond}, '_', '\\_'), ' ', currSummary{iCond, 3}{:}, '  (n = ', num2str(currSummary{iCond, 1}), ')'];
    end
   
    % Plot figures
    dffCurrConds = combinedDffAvg{eventInd}(:,:,:, currConds);
    range = calc_range(dffCurrConds, rangeScalar, rangeType);
    [~, ~] = plot_heatmaps(dffCurrConds, infoStruct, range, plotTitles(currConds), sigma, 'fileName', fileName, 'makeVid', makeVid, 'saveDir', saveDir);
     
clear odorEventName eventInd currSummary currCondNames currConds sigma rangeType rangeScalar makeVid saveDir fileName plotTitles dffCurrConds range
catch foldME; rethrow(foldME); end  
%% =================================================================================================
%   BEHAVIOR HEATMAPS                                   
%%%=================================================================================================
    %% PLOT OVERALL MEAN dF/F ACROSS BEHAVIORAL STATES
try
    behaviorNames = {'Locomotion', 'IsoMove', 'AnyMove'};
    actionNum = [1];

    smoothingSigma = [0.6]; 
    rangeType = 'max';
    rangeScalar = 0.5;
    makeVid = 1;
    saveDir = [];
    fileName = ['All_frame_dFF_', behaviorNames{actionNum}, '_Heatmaps_', regexprep(expDate, {'_', 'exp'}, {'', '-'})];
    titleStr = {['dF/F - ', behaviorNames{actionNum}, ' vs. Quiescence']};
    
    % Load dF/F data
    load(fullfile(expDir, ['actionDff_', behaviorNames{actionNum}, '.mat'])) % --> 'actionDff'

    % Calculate absolute max dF/F value across all planes and action states
    range = calc_range(actionDff, rangeScalar, rangeType);

    % Plot figures
    [f, ~] = plot_heatmaps(actionDff, infoStruct, range, titleStr, smoothingSigma, 'fileName', fileName, 'makeVid', makeVid, ...
                           'saveDir', saveDir);
clear smoothingSigma rangeType rangeScalar makeVid saveDir fileName titleStr range   
catch foldME; rethrow(foldME); end
    %% PLOT INTERACTION HEATMAPS FOR SOME TRIAL CONDITIONS
try
    % Show summary again
    eventInd = contains(primaryEventNames, 'locomotion');
    currSummary = allCondSummaries{eventInd};
    disp(currSummary)
                 
    currCondNames = repmat(allCondNames{eventInd}, 2, 1);
    
    % Calculate absolute max dF/F value across all planes and stim types
    currConds = [2 4];
    sigma = [0.2]; 
    rangeType = 'Max';
    rangeScalar = 0.6;
    makeVid = 0;
    saveDir = [];
    fileName = ['Locomotion_Response_Heatmaps_', regexprep(expDate, {'_', 'exp'}, {'', '-'})];

    plotTitles = [];
    for iCond = currConds
        plotTitles{iCond} = [regexprep(currCondNames{iCond}, '_', '\\_'), ' ', currSummary{iCond, 3}{:}, '  (n = ', num2str(currSummary{iCond, 1}), ')'];
    end
   
    % Plot figures
    dffCurrConds = combinedDffAvg{eventInd}(:,:,:, currConds);
    range = calc_range(dffCurrConds, rangeScalar, rangeType);
    [~, ~] = plot_heatmaps(dffCurrConds, infoStruct, range, plotTitles(currConds), sigma, 'fileName', fileName, 'makeVid', makeVid, 'saveDir', saveDir); 

    clear eventInd currSummary currCondNames currConds sigma rangeType rangeScalar makeVid saveDir fileName plotTitles range

catch foldME; rethrow(foldME); end

%% =================================================================================================
%           ROI-BASED ANALYSES                                   
%%%=================================================================================================
    %% PLOT AND SAVE NEW ROIs
try
ROIselectionGui();
catch foldME; rethrow(foldME); end
    %% LOAD ROI DATA AND PLOT ROIs
try
activeROIs = [];
expDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)];
metaDataFileName = 'ROI_metadata.mat';
dffDataFileName = 'ROI_Data_Avg.mat';


% ROInames = ["L-SLP", "R-SLP", "L-ANT", "R-ANT", "Control"];
% ROInames = ["L-SLP", "R-SLP", "L-ANT", "R-ANT","Bi-ANT", "Control"];
% ROInames = ["SLP", "ANT", "PPL1", "Control"];
ROInames = ["SLP", "ANT", "LH", "PPM2", "Control"];
% ROInames = ["SLP", "ANT", "LH", "PPM2", "VLP", "SMP", "Control"];
% ROInames = ["L-SLP", "R-SLP", "Bi-ANT", "Control"];
% metaDataFileName = 'ROI_metadata_SLP_comparison.mat';
% dffDataFileName = 'ROI_Data_Avg_SLP_comparison.mat';

activeROIs = [];

% Load metadata
load(fullfile(expDir, metaDataFileName)); % --> ROImetadata(.mask, .xi, .yi, .plane, .color, .refImg)
if ~isempty(activeROIs)
    ROImetadata = ROImetadata(activeROIs);
end
infoStruct.ROImetadata = ROImetadata;
infoStruct.nROIs = numel(ROImetadata); nROIs = infoStruct.nROIs;

% Load imaging and dF/F data
load(fullfile(expDir, dffDataFileName)); % --> 'ROIDataAvg', 'ROIDffAvg', 'ROIDataBaseSub' ([volume, trial, ROI])
if ~isempty(activeROIs)
    ROIDataAvg = ROIDataAvg(:,:, activeROIs);
    ROIDffAvg = ROIDffAvg(:,:, activeROIs);
    if max(activeROIs) <= size(ROIDataBaseSub, 3)
        ROIDataBaseSub = ROIDataBaseSub(:,:, activeROIs);
    else
        currROIs = activeROIs(activeROIs ~= max(activeROIs));
        ROIDataBaseSub = ROIDataBaseSub(:,:,currROIs);
    end
end

disp('ROI data loaded')

% Create ROI names if necessary
if ~isempty(activeROIs) && exist('ROInames', 'var')
    ROInames = ROInames(activeROIs);
end
if ~exist('ROInames', 'var') || numel(ROInames) ~= nROIs
    ROInames = join([repmat("ROI", 1, nROIs); string(1:nROIs)], 1);
end
for iROI = 1:numel(ROImetadata)
    try
        ROImetadata{iROI}(1).name = ROInames{iROI};
    catch
        ROImetadata{iROI}(1).name = '';
    end
end

% Plot ROIs 
close all;
plot_ROIs(ROImetadata, 'IntensityRange', [0 MAX_INTENSITY]);

% Plot mean value in each ROI for each trial with and without subtracting last ROI
volAvgROIData = squeeze(mean(ROIDataAvg, 1)); % --> [trial, ROI]
f = figure(nROIs + 1);clf;
f.Position(3) = 2 * f.Position(4);
subaxis(1,2,1);hold on
legendStr = [];
for iROI = 1:nROIs
    currROIData = volAvgROIData(:, iROI);
    zeroedData = currROIData - min(currROIData(:));
    legendStr{iROI} = ROInames{iROI};
    plot(zeroedData);
end
legend(legendStr);
subaxis(1,2,2);hold on
legendStr = [];
for iROI = 1:nROIs
    currROIData = volAvgROIData(:, iROI);
    lastROIData = volAvgROIData(:, end);
    zeroedData = currROIData - min(currROIData(:));
    lastROIZeroed = lastROIData - min(lastROIData(:));
    zeroedData = zeroedData - lastROIZeroed;
    legendStr{iROI} = ROInames{iROI};
    plot(zeroedData);
end
legend(legendStr);


clear metaDataFileName dffDataFileName 
catch foldME; rethrow(foldME); end
    %% PLOT 2D HEATMAPS OF ROI FLUORESCENCE THROUGHOUT EXPERIMENT 
try
saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), ...
        '\Analysis'], 'Select a save directory');
sepBlockStims = 0; sepGroupStims = 0; clear minMax;
fontSize = 14;

useDff = 0;
subtractBaseline = 0;
colorbarOn = 1;
% minMax = [-0.12 1.5; -0.12 1; -0.12 0.8; -0.12 2.5; 0 1]; % [ROI, min-max]
% minMax = [110 270; 300 750; 150 1000; 75 425; 0 600; 50 700]; % [ROI, min-max]
% % % 

% trialGroups = [];
% plotTitleSuffix = 'All Trials';
% fileNameSuffix = '_AllTrials';
% 
% % % % % % % % % % % % plotTitleSuffix = ['sid\_', num2str(sid)]

% % GROUP BY STIM TYPE
% trialGroups = stimTrialGroups; 
% plotTitleSuffix = make_plotTitleSuffix(stimNames); %
% fileNameSuffix = make_fileNameSuffix(stimGroupNames);
% sepGroupStims = 1;

% GROUP CHRONOLOGICALLY BY BLOCKS
trialGroups = zeros(1, nTrials);
for iBound = 1:numel(groupBounds)-1
   trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
end
trialGroups(groupBounds(end)+1:end) = iBound + 1;
trialGroups = trialGroups .* goodTrials;
plotTitleSuffix = '';
fileNameSuffix = '_Blocks_Separated';
if ~isempty(blockShading)
    sepBlockStims = 1;
end


%---------------------------------------------------------------------------------------------------
try
% Create base file name and select data
if useDff
    fileNamePrefix = ['ROI_2D_dFF_Summary_', regexprep(expDate, {'_', 'exp'}, {'', '-'})];
    flData = ROIDffAvg;
    baselineData = ROIDffAvg(:,:,end);
else
    fileNamePrefix = ['ROI_2D_Fluorescence_Summary_', regexprep(expDate, {'_', 'exp'}, {'', '-'})];
    flData = ROIDataAvg;
    baselineData = ROIDataAvg(:,:,end);
end

% Subtract baseline if necessary
if subtractBaseline
    currNumROIs = nROIs - 1;
    fileNamePrefix = [fileNamePrefix, '_BaselineSub'];
    flData = flData - repmat(baselineData, 1, 1, size(flData, 3));
else
    currNumROIs = nROIs;
end

for iROI = 1:currNumROIs
    
    % Update plot title
    if useDff
        plotTitlePrefix = [regexprep(expDate, {'_', 'exp'}, {'', '-'}), ' — ', ROInames{iROI}, ...
            '  dF/F '];
    else
        plotTitlePrefix = [regexprep(expDate, {'_', 'exp'}, {'', '-'}), ' — ', ROInames{iROI}, ...
            '  raw fluorescence'];
    end
    if subtractBaseline
        plotTitlePrefix = [plotTitlePrefix, '  (baseline-sub)'];
    end
    
    
    % Create figure
    f = figure(iROI); clf
    f.Color = [1 1 1];
    f.Position = [-1050 10 1020 950];
    ax = axes();
    
    % Create colormap
    currFlData = flData(:,:,iROI);
    if exist('minMax', 'var')
        if ~isnan(minMax(iROI, 1))
            currFlData(1,1) = minMax(iROI, 1);
            currFlData(currFlData < minMax(iROI, 1)) = minMax(iROI, 1);
        end
        if ~isnan(minMax(iROI, 2))
            currFlData(end,end) = minMax(iROI, 2);
            currFlData(currFlData > minMax(iROI, 2)) = minMax(iROI, 2);
        end
    end
    cm = parula(numel(unique(currFlData(:,goodTrials))));
    currFlData(:, ~goodTrials) = max(as_vector(currFlData(:,goodTrials))) * 1.01;
    if sum(~goodTrials)
        cm = [0 0 0; cm(2:end-1, :); 1 1 1];
    end
    
    % Plot fluorescence heatmap
    plot_2D_summary(currFlData', volumeRate, ...
        'plotAxes', ax, ...
        'trialGroups', trialGroups, ...
        'titleStr', [plotTitlePrefix, ' — ', plotTitleSuffix], ...
        'colormap', cm ...
        );
    if colorbarOn
        colorbar;
    end
    if sepBlockStims
        ax.YTick = [];
    end
    ax.FontSize = fontSize;
    
% Plot stim times
hold on
if sepBlockStims
    
    % Calculate y-axis ranges for each block
    blockStarts = 1;
    blockEnds = [];
    for iBlock = 1:numel(unique(trialGroups(trialGroups > 0)))-1
        blockTrials = sum(trialGroups == iBlock);
        blockStarts(iBlock + 1) = blockStarts(iBlock) + 4 + blockTrials;
        blockEnds(iBlock) = blockStarts(iBlock) + blockTrials - 1;
    end
    blockEnds(end + 1) = numel(trialGroups) + (4 * numel(unique(trialGroups(trialGroups > 0))));
    blockRanges = [blockStarts - 1; blockEnds + 1]';
    
    % Plot stim start and end for each block
    for iBlock = 1:numel(unique(trialGroups(trialGroups > 0)))
        currBlockShading = stimEpochs(blockShading{iBlock}, :);
        nStimEpochs = size(currBlockShading, 1);
        epochInds = blockShading{iBlock};
        currShadeColors = stimShadingColors(blockShading{iBlock});
        draw_stim_lines(currBlockShading, currShadeColors, 'plotAxes', ax, 'yLims', ...
            blockRanges(iBlock, :), 'frameRate', volumeRate);
    end
elseif sepGroupStims
    
    % Calculate y-axis ranges for each stim type
    blockStarts = 1;
    blockEnds = [];
    for iBlock = 1:numel(unique(trialGroups(trialGroups > 0)))-1
        blockTrials = sum(trialGroups == iBlock);
        blockStarts(iBlock + 1) = blockStarts(iBlock) + 4 + blockTrials;
        blockEnds(iBlock) = blockStarts(iBlock) + blockTrials - 1;
    end
    blockEnds(end + 1) = numel(trialGroups) + (4 * numel(unique(trialGroups(trialGroups > 0))));
    blockRanges = [blockStarts - 1; blockEnds + 1]';
    
    % Plot stim start and end for each block
    for iBlock = 1:numel(unique(trialGroups(trialGroups > 0)))
        currBlockShading = stimEpochs(stimGroupShading{iBlock}, :);
        nStimEpochs = size(currBlockShading, 1);
        epochInds = blockShading{iBlock};
        currShadeColors = stimShadingColors(stimGroupShading{iBlock});
        draw_stim_lines(currBlockShading, currShadeColors, ...
            'plotAxes', ax, 'yLims', blockRanges(iBlock, :), 'frameRate', volumeRate);
    end
else
    % Draw one set of lines for the entire plot
    draw_stim_lines(stimEpochs, stimShadingColors, 'plotAxes', ax, 'frameRate', volumeRate);
end

% Save figure
if saveDir
    fileName = [fileNamePrefix, '_', strrep(ROInames{iROI}, ' ', '_'), fileNameSuffix];
    save_figure(f, saveDir, fileName);
end

end%iROI
catch foldME; rethrow(foldME); end    
catch foldME; rethrow(foldME); end
    %% EXTRACT EVENT-TRIGGERED dF/F WITHIN ROIs
try
ROIEventDff = [];
for iType = 1:nEventTypes
    
    primaryFiltName = primaryEventNames{iType};
    analysisWindow = analysisWindows(iType, :);
    baselineDur = analysisWindow(1);
    respDur = analysisWindow(2);
    
    currCondSum = allCondSummaries{iType};
    nConds = numel(currCondSum.CondName);
        
    trialNums = 1:nTrials;
    goodTrialNums = trialNums(goodTrials);
    goodEvents = ismember(eventLists{iType}(:,3), goodTrialNums);
    eventList = eventLists{iType}(goodEvents,:);
    
    for iCond = 1:nConds
        
        disp(['Extracting event data for ', primaryFiltName, ' cond #', num2str(iCond), ' of ', num2str(nConds), '...'])
                
        % Get ROI data for event onsets
        offsetAlign = strcmp(currCondSum.Align{iCond}, 'offset');
        [baselineData, respData] = extract_event_volumes(eventList, combFilterVecs{iType}(goodEvents,iCond), baselineDur, respDur, infoStruct, ...
            permute(ROIDffAvg, [3 1 2]), 'offsetAlign', offsetAlign);   % --> [ROI, volume, event]
        baselineData = permute(baselineData, [2 3 1]);                  % --> [volume, event, ROI]
        respData = permute(respData, [2 3 1]);                          % --> [volume, event, ROI]
        combDff = cat(1, baselineData, respData);                       % --> [volume, event, ROI]
        ROIEventDff{iType}{iCond} = combDff;                            % --> {eventType}{cond}[volume, event, ROI]
        
    end%iCond
    clear respData baselineData combDff
    disp('ROI event data extracted');
end% iType
clear primaryFiltName analysisWindow eventList baselineDur respDur currMask currPlan nEvents baselineVols
catch foldME; rethrow(foldME); end
            %% PLOT EVENT-ALIGNED dF/F WITHIN ROIs
try
% Show summary again
eventName = 'OdorA';
fileNamePrefix = 'Locomotion_responses_';
shadeDur = 0;
eventInd = contains(primaryEventNames, eventName);
currSummary = allCondSummaries{eventInd};
disp(currSummary)

currConds = [2 4];
currCondNames = allCondNames{eventInd}(currConds);

currDffData = ROIEventDff{eventInd}(currConds); % --> {cond}[volume, event, ROI]

saveDir = 0;
saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), ...
        '\Analysis'], 'Select a save directory');


ROIlist = 1:size(currDffData{1}, 3)-1;
% ROIlist = [1 2 3];

ax = []; yLims = [];
for iROI = ROIlist
    
    % Create figure
    f = figure(iROI); clf
    f.Position = [-1850 100 1200 800];
    f.Color = [1 1 1];
    
    % Plot reference image and ROI outline
    plotPlaneNum = 1;
    currPlane = ROImetadata{iROI}(plotPlaneNum).plane;
    xData = ROImetadata{iROI}(plotPlaneNum).xi;
    yData = ROImetadata{iROI}(plotPlaneNum).yi;
    nPlots = numel(currDffData) + 1;
    if nPlots == 3
        plotPos = [2 2]; % Because [1 3] looks bad
    else
        plotPos = numSubplots(nPlots);
    end
    subaxis(plotPos(1), plotPos(2), 1,'ML', 0.05, 'MR', 0.02, 'MT', 0.05, 'MB', 0.08, 'SV', 0.1, ...
            'SH', 0.05)
    hold on
    imshow(infoStruct.refImg{currPlane}, [0 MAX_INTENSITY]);
    plot(xData, yData, 'Color', 'g');
    title(ROInames{iROI});
%     title(['Plane #', num2str(ROImetadata{iROI}(1).plane)])
    
    % Plot dF/F for each condition
    alignStr = currSummary.Align(currConds);
    for iPlot = 1:(nPlots - 1) 
        if strcmp(alignStr{iPlot}, 'onset')
            eventShading = [0, shadeDur];
        else
            eventShading = [-shadeDur, 0];
        end
        plotDff = currDffData{iPlot}(:,:,iROI);
        if nPlots == 3
            subaxis(plotPos(1), plotPos(2), iPlot + 2); % If there's only two plots they look better side-by-side
        elseif nPlots == 5
            subaxis( ((mod(iPlot,2) + iPlot) / 2) + iPlot ); % If there's four they should be in a square
        else
            subaxis(plotPos(1), plotPos(2), iPlot + 1)
        end
        ax{iPlot} = gca;
        volOffset = round(analysisWindows(eventInd, 1) * volumeRate * -1);
        plot_ROI_data(ax{iPlot}, plotDff, 'VolumeRate', volumeRate, 'VolOffset', volOffset, ...
                'OutlierSD', 5); 
        title(regexprep([expDate, ' — ', ROInames{iROI}, '  —  ', currCondNames{iPlot}, ' - ', ...
                alignStr{iPlot}], '_', '\\_'))
        yLims(iPlot,:) = ylim(ax{iPlot});
    end
    
    % Scale y-axes to match
    yMin = min(yLims(:));
    yMax = max(yLims(:));
    for iPlot = 1:(nPlots - 1)
       ylim(ax{iPlot}, [yMin yMax]); 
    end
    
    % Save figure ----------------------------------------------------------------------------------
    if saveDir
        fileName = [fileNamePrefix, strrep(ROInames{iROI}, ' ', '_'), '_', expDate];
        export_fig(fullfile(saveDir, fileName), '-png', f);
        if ~isdir(fullfile(saveDir, 'figFiles'))
            mkdir(fullfile(saveDir, 'figFiles'))
        end
        savefig(f, fullfile(saveDir, 'figFiles', fileName));
    end
end
catch foldME; rethrow(foldME); end
    %% PLOT MEAN ROI dF/F THROUGHOUT TRIAL IN SEPARATE PLOTS FOR EACH TRIAL GROUP
try
    saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', ...
            num2str(sid), '\Analysis'], 'Select a save directory');
    trialGroups = goodTrials;
    cm = []; sepBlockStims = 0; fileNameSuffix = '';
    figPos = [-1450 50 900 700]; %[ -1450 50 800 900];% 
    
    useDff = 0;
    subtractBaseline = 1;
    
    smWin = 3;
    singleTrials = 0;
    singleTrialAlpha = 0.5;
    stdDevShading = 1;
    outlierSD = 4;

%     % ALL TRIALS
%     trialGroups = [goodTrials];
%     fileNameSuffix = '_AllTrials';
%     plotTitleSuffix = '  —  All Trials';
%     %     trialGroups(31:end) = 0;
    
%     % GROUP CHRONOLOGICALLY BY BLOCKS
%     trialGroups = zeros(1, nTrials);
%     for iBound = 1:numel(groupBounds)-1
%         trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
%     end
%     trialGroups(groupBounds(end)+1:end) = iBound + 1;
%     trialGroups = trialGroups .* goodTrials;
%     plotTitleSuffix = '';
%     fileNameSuffix = '_Blocks_Separated';
%     if ~isempty(blockShading)
%         sepBlockStims = 1;
%     end
    %     plotTitleSuffix = ['  —  sid\_', num2str(sid)]

    % GROUP BY STIM TYPE
    trialGroups = stimTrialGroups .* goodTrials;
    fileNameSuffix = '_StimTypeComparison';
    plotTitleSuffix = make_plotTitleSuffix(stimNames);
    
% --------------------------------------------------------------------------------------------------
try
% Create base file name and select data
if useDff
    fileNamePrefix = ['Trial_Avg_Dff_', regexprep(expDate, {'_', 'exp'}, {'', '-'})];
    flData = ROIDffAvg;
    baselineData = ROIDffAvg(:,:,end);
else
    fileNamePrefix = ['Trial_Avg_RawF_', regexprep(expDate, {'_', 'exp'}, {'', '-'})];
    flData = ROIDataAvg;
    baselineData = ROIDataAvg(:,:,end);
end

% Subtract baseline if necessary
if subtractBaseline
    currNumROIs = nROIs - 1;
    fileNamePrefix = [fileNamePrefix, '_BaselineSub'];
    flData = flData - repmat(baselineData, 1, 1, size(flData, 3));
    flData(:,:,end) = [];
else
    currNumROIs = nROIs;
end

    ROIlist = 1:currNumROIs;
%     ROIlist = [1 2 3];
    yL = [];
    for iROI = ROIlist       

        nGroups = numel(unique(trialGroups(trialGroups > 0)));
        nRows = nGroups + 1;
        
        % Update plot title and axis label
        if useDff
            plotTitlePrefix = [regexprep(expDate, {'_', 'exp'}, {'', '-'}), ...
                    '    ', ROInames{iROI}];
            yAxisLabel = 'dF/F';
        else
            plotTitlePrefix = [regexprep(expDate, {'_', 'exp'}, {'', '-'}), ...
                     '    ', ROInames{iROI}];
            yAxisLabel = 'Raw F';
        end
        if subtractBaseline
            plotTitlePrefix = [plotTitlePrefix, '  (baseline-sub)'];
        end

        % Create figure
        f = figure(iROI); clf
        f.Position = figPos;
        f.Color = [1 1 1];
        
        % Plot reference image and ROI outline
        currPlane = ROImetadata{iROI}(1).plane;
        xData = ROImetadata{iROI}(1).xi;
        yData = ROImetadata{iROI}(1).yi;
        subaxis(nRows, 2, [1 2], 'ML', 0.06, 'MR', 0.015, 'MT', 0.05, 'MB', 0.08, 'SH', 0)
        hold on
        imshow(infoStruct.refImg{currPlane}, [0 infoStruct.MAX_INTENSITY]);
        plot(xData, yData, 'Color', 'g');
        title(ROInames{iROI})
        
        clear ax
        for iGroup = 1:nGroups
            currFlData = flData(:, trialGroups == iGroup, iROI); % --> [volume, trial]
            subaxis(nRows, 2, 1, iGroup+1, 2,1)  %( (iPlot*3 + 1):(iPlot*3 + 3) ))                                                     
            ax(iGroup) = gca;
            ax(iGroup) = plot_ROI_data(ax(iGroup), currFlData, ...
                                            'SingleTrials', singleTrials, ...
                                            'SingleTrialAlpha', singleTrialAlpha, ...
                                            'StdDevShading', stdDevShading, ...
                                            'OutlierSD', outlierSD, ...
                                            'VolumeRate', volumeRate, ...
                                            'YAxisLabel', yAxisLabel, ...
                                            'SmoothWinSize', smWin);
                                        
            yL{iGroup} = ylim(ax(iGroup));
            if iGroup ~= nGroups
                xlabel([]);
            end
            hold on
            if sepBlockStims
                % Plot stim start(s) and end(s) for the current block
                currBlockShading = blockShading{iGroup};
                currStimEpochs = stimEpochs(currBlockShading, :);
                currStimShadingColors = stimShadingColors(currBlockShading);
            else
                % Just draw one set of lines for each stim epoch
                currStimEpochs = stimEpochs;
                currStimShadingColors = stimShadingColors;
            end
            for iEpoch = 1:size(currStimEpochs, 1)
                plot_stim_shading(currStimEpochs(iEpoch, :), 'Axes', ax(iGroup), 'Color', ...
                    rgb(currStimShadingColors{iEpoch}));
            end
        end
        ax(1).Title.String = [plotTitlePrefix, plotTitleSuffix];

        % Scale y-axes to match
        yMin = min([yL{:}]);
        yMax = max([yL{:}]);
        for iGroup = 1:nGroups
            ylim(ax(iGroup), [yMin yMax]);
        end
        
        % Save figure ------------------------------------------------------------------------------
        if saveDir
            fileName = [fileNamePrefix, fileNameSuffix, '_', strrep(ROInames{iROI}, ' ', '_')];
            save_figure(f, saveDir, fileName);
        end
    end
catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end 
    %% PLOT MEAN ROI dF/F THROUGHOUT TRIAL OVERLAYING MULTIPLE TRIAL GROUPS
try
    saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', ...
            num2str(sid), '\Analysis'], 'Select a save directory');
    trialGroups = goodTrials;
    cm = []; fileNameSuffix = '';
    figPos = [ -1450 50 800 900];%[-1450 50 1000 900] 
    
    useDff = 0;
    subtractBaseline = 1;
    
    smWin = 3;
    singleTrials = 0;
    singleTrialAlpha = 0.3;
    stdDevShading = 0;
    outlierSD = 5;

%     % GROUP BY STIM TYPE
%     trialGroups = stimTrialGroups .* goodTrials;
%     fileNameSuffix = '_StimTypeComparison';
%     binNames = stimGroupNames;

%     % GROUP CHRONOLOGICALLY BY BLOCKS
%     trialGroups = zeros(1, nTrials);
%     for iBound = 1:numel(groupBounds)-1
%         trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
%     end
%     trialGroups(groupBounds(end)+1:end) = iBound + 1;
%     trialGroups = trialGroups .* goodTrials;
%     plotTitleSuffix = '';
%     fileNameSuffix = '_Blocks_Separated';
%     if ~isempty(blockShading)
%         sepBlockStims = 1;
%     end
%     binNames = blockNames;
    
%         % Skip first block
%         trialGroups = trialGroups - 1;
%         trialGroups(1:20) = 0;
%         binNames = binNames(2:4);

    % GROUP BY BLOCK TYPE
    groupNames = blockNames(1:2);
    trialGroups = zeros(1, nTrials);
    trialGroups = trialGroups .* goodTrials;
    for iBound = 1:numel(groupBounds)-1
        trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
    end
    trialGroups(logical(mod(trialGroups, 2))) = 1;
    trialGroups(~logical(mod(trialGroups, 2))) = 2;
    fileNameSuffix = '_Block_Type_Groups';

%     % GROUP BY EARLY/LATE
%     fileNameSuffix = '_EarlyVsLateTrials';
%     plotTitleSuffix = '  —  Early vs. late trials';
%     binNames = [];
%     binBounds = [1, 20];
%     trialGroups = zeros(1, nTrials);
%     for iBound = 1:numel(binBounds)-1
%         trialGroups(binBounds(iBound):binBounds(iBound + 1)) = iBound;
%         binNames{iBound} = ['Trials ', num2str(binBounds(iBound)), '-', num2str(binBounds(iBound + 1))];
%     end
%     trialGroups(binBounds(end):end) = numel(binBounds);
%     trialGroups = trialGroups .* goodTrials;
%     binNames{end + 1} = ['Trials ', num2str(binBounds(end)), '-', num2str(nTrials)];

    % --------------------------------------------------------------------------------------------------
try
    % Create base file name and select data
    if useDff
        fileNamePrefix = ['Trial_Avg_Overlay_Dff_', regexprep(expDate, {'_', 'exp'}, {'', '-'})];
        flData = ROIDffAvg;
        baselineData = ROIDffAvg(:,:,end);
        yAxisLabel = 'dF/F';
    else
        fileNamePrefix = ['Trial_Avg_Overlay_RawF_', regexprep(expDate, {'_', 'exp'}, {'', '-'})];
        flData = ROIDataAvg;
        baselineData = ROIDataAvg(:,:,end);
        yAxisLabel = 'Raw F';
    end

    % Subtract baseline if necessary
    if subtractBaseline
        currNumROIs = nROIs - 1;
        fileNamePrefix = [fileNamePrefix, '_BaselineSub'];
        flData = flData - repmat(baselineData, 1, 1, size(flData, 3));
        flData(:,:,end) = [];
    else
        currNumROIs = nROIs;
    end

    ROIlist = 1:currNumROIs;
%     ROIlist = [1 2 3];
    for iROI = ROIlist       

        % Create figure
        f = figure(iROI); clf
        f.Position = figPos;
        f.Color = [1 1 1];

        % Plot reference image and ROI outline
        currPlane = ROImetadata{iROI}(1).plane;
        xData = ROImetadata{iROI}(1).xi;
        yData = ROImetadata{iROI}(1).yi;
        subaxis(2, 3, [1 2], 'ML', 0.06, 'MR', 0.015, 'MT', 0.05, 'MB', 0.08, 'SH', 0)
        hold on
        imshow(infoStruct.refImg{currPlane}, [0 infoStruct.MAX_INTENSITY]);
        plot(xData, yData, 'Color', 'g');
        title(ROInames{iROI})
        
        % Create plot
        currFlData = flData(:,:, iROI); % --> [volume, trial]
        subaxis(2, 3, 4:6)
        ax = gca;
        ax = plot_ROI_data(ax, currFlData, ...
                                'TrialGroups', trialGroups,   ...
                                'SingleTrials', singleTrials, ...
                                'SingleTrialAlpha', singleTrialAlpha, ...
                                'OutlierSD', outlierSD, ...
                                'Legend', groupNames, ...
                                'VolumeRate', volumeRate, ...
                                'YAxisLabel', yAxisLabel, ...
                                'StdDevShading', stdDevShading, ...
                                'SmoothWinSize', smWin);
        title([regexprep(expDate, {'_', 'exp'}, {'', '-'}), ...
                '   ', ROInames{iROI}, plotTitleSuffix]);
        
        % Plot stim shading epoch(s)
        for iEpoch = 1:size(stimEpochs, 1)
            plot_stim_shading(stimEpochs(iEpoch, :), 'Axes', ax, 'Color', ...
                rgb(stimShadingColors{iEpoch}));
        end

        % Save figure -------------------------------------------------------------------------------
        if saveDir
            fileName = [fileNamePrefix, fileNameSuffix, '_', strrep(ROInames{iROI}, ' ', '_')];
            save_figure(f, saveDir, fileName);
        end
    end
catch foldME; rethrow(foldME); end    
catch foldME; rethrow(foldME); end
    %% PLOT MEAN ROI dF/F THROUGHOUT TRIAL COLOR CODED BY ANNOTATION DATA
try
 saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', ...
        num2str(sid), '\Analysis'], 'Select a save directory');
trialGroups = goodTrials;
sepBlockStims = 0; fileNameSuffix = '';
edgeColorMode = 'flat';
cm =[rgb('RoyalBlue') * 0.6; ...
    rgb('Orange'); ...
    rgb('Crimson'); ...
    rgb('Cyan') * 0.85];
figPos = [-1450 50 1000 900]; %[ -1450 50 800 900];%

useDff = 0;
subtractBaseline = 1;
annotValues = [3 0]; 
cm = cm(1:max(annotValues) + 1, :);
annotTypeInd = 4;

smWin = 3;
singleTrials = 1;
singleTrialAlpha = 0.4;
stdDevShading = 0;
outlierSD = 4;

behavData = annotationTypes{annotTypeInd}.volAnnotArr';
behavData(~ismember(behavData, annotValues)) = nan;

% ALL TRIALS
trialGroups = goodTrials+0;
fileNameSuffix = '_AllTrials';
plotTitleSuffix = '  —  All Trials';
% 
% plotTitleSuffix = ['  —  sid\_', num2str(sid)]
    
% % GROUP BY STIM TYPE
% trialGroups = stimTrialGroups .* goodTrials;
% fileNameSuffix = '_StimTypeComparison';
% plotTitleSuffix = make_plotTitleSuffix(stimNames);

% --------------------------------------------------------------------------------------------------
try
% Create base file name and select data
if useDff
    fileNamePrefix = ['Behavior_Annot_Coded_Avg_Dff_', regexprep(expDate, {'_', 'exp'}, {'', '-'})];
    flData = ROIDffAvg;
    baselineData = ROIDffAvg(:,:,end);
else
    fileNamePrefix = ['Behavior_Annot_Coded_Avg_RawF_', regexprep(expDate, {'_', 'exp'}, {'', '-'})];
    flData = ROIDataAvg;
    baselineData = ROIDataAvg(:,:,end);
end

% Subtract baseline if necessary
if subtractBaseline
    currNumROIs = nROIs - 1;
    fileNamePrefix = [fileNamePrefix, '_BaselineSub'];
    flData = flData - repmat(baselineData, 1, 1, size(flData, 3));
    flData(:,:,end) = [];
else
    currNumROIs = nROIs;
end

ROIlist = 1:currNumROIs;
%     ROIlist = [1 2 3];
yL = []; %f = [];
for iROI = ROIlist       

    nGroups = numel(unique(trialGroups(trialGroups > 0)));
    nRows = nGroups + 1;

    % Update plot title
    if useDff
        plotTitlePrefix = [regexprep(expDate, {'_', 'exp'}, {'', '-'}), ...
                '    ', ROInames{iROI}];
        yAxisLabel = 'dF/F';
    else
        plotTitlePrefix = [regexprep(expDate, {'_', 'exp'}, {'', '-'}), ...
                 '    ', ROInames{iROI}];
        yAxisLabel = 'Raw F';
    end
    if subtractBaseline
        plotTitlePrefix = [plotTitlePrefix, '  (baseline-sub)'];
    end

    % Create figure
    f = figure(iROI); clf
    f.Position = figPos;
    f.Color = [1 1 1];

    % Plot reference image and ROI outline
    currPlane = ROImetadata{iROI}(1).plane;
    xData = ROImetadata{iROI}(1).xi;
    yData = ROImetadata{iROI}(1).yi;
    subaxis(nRows, 2, [1 2], 'ML', 0.06, 'MR', 0.015, 'MT', 0.05, 'MB', 0.08, 'SH', 0)
    hold on
    imshow(infoStruct.refImg{currPlane}, [0 infoStruct.MAX_INTENSITY]);
    plot(xData, yData, 'Color', 'g');
    title(ROInames{iROI})
    
    clear ax
    for iGroup = 1:nGroups
        currFlData = flData(:, trialGroups == iGroup, iROI); % --> [volume, trial]
        currBehavData = behavData(:, trialGroups == iGroup);
        subaxis(nRows, 2, 1, iGroup+1, 2,1)  %( (iPlot*3 + 1):(iPlot*3 + 3) ))                                                     
        ax(iGroup) = gca;
        ax(iGroup) = plot_behavior_coded_ROI_data(ax(iGroup), currFlData, currBehavData, ...
                   'SingleTrials', singleTrials, ...
                   'SingleTrialAlpha', singleTrialAlpha, ...
                   'StdDevShading', stdDevShading, ...
                   'OutlierSD', outlierSD, ...
                   'VolumeRate', volumeRate, ...
                   'EdgeColorMode', edgeColorMode, ...
                   'Colormap', cm, ...
                   'SmoothWinSize', smWin);
                   
        yL{iGroup} = ylim(ax(iGroup));
        if ~useDff
            ylabel('Raw F');

        end
        if iGroup ~= nGroups 
            xlabel([]);
        end
        hold on
        if sepBlockStims
            % Plot stim start(s) and end(s) for the current block
            currBlockShading = blockShading{iGroup};
            currStimEpochs = stimEpochs(currBlockShading, :);
            currStimShadingColors = stimShadingColors(currBlockShading);
        else
            % Just draw one set of lines for each stim epoch
            currStimEpochs = stimEpochs;
            currStimShadingColors = stimShadingColors;
        end
        for iEpoch = 1:size(currStimEpochs, 1)
            plot_stim_shading(currStimEpochs(iEpoch, :), 'Axes', ax(iGroup), 'Color', ...
                rgb(currStimShadingColors{iEpoch}));
        end
    end
    ax(1).Title.String = [plotTitlePrefix, plotTitleSuffix];
    
    % Exclude first and last volumes to improve appearance
    xL = xlim;
    xlim([xL(1) + 1/volumeRate, xL(2) - 1/volumeRate])
    
    % Scale y-axes to match
    yMin = min([yL{:}]);
    yMax = max([yL{:}]);
    for iGroup = 1:nGroups
        ylim(ax(iGroup), [yMin yMax]);
    end

    % Save figure ------------------------------------------------------------------------------
    if saveDir
        fileName = [fileNamePrefix, fileNameSuffix, '_', strrep(ROInames{iROI}, ' ', '_')];
        save_figure(f, saveDir, fileName);
    end
end
catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end
    %% PLOT MEAN ROI dF/F THROUGHOUT TRIAL COLOR CODED BY FICTRAC DATA
try
 saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', ...
        num2str(sid), '\Analysis'], 'Select a save directory');
trialGroups = goodTrials;
cm = []; groupShading = []; fileNameSuffix = ''; CLimCap = 1;
figPos = [ -1450 50 800 900]; %[-1450 50 1000 900]; %
fontSize = 8;

useDff = 0;
subtractBaseline = 1;

smWin = 5;
singleTrials = 0;
singleTrialAlpha = 0.4;
stdDevShading = 0;
outlierSD = 5;

% Color by FicTrac data
edgeColorMode = 'interp';
% cm = 'bluewhitered';
cm = 'parula';
ftVar = 'moveSpeed';%'yawVel';%'yawSpeed';%   
CLimCap = 1;

% % ALL TRIALS
% trialGroups = goodTrials;
% fileNameSuffix = '_AllTrials';
% plotTitleSuffix = '  —  All Trials';

% plotTitleSuffix = ['  —  sid\_', num2str(sid)]
% 
% % GROUP BY STIM TYPE
% trialGroups = stimTrialGroups .* goodTrials;
% fileNameSuffix = '_StimTypeComparison';
% plotTitleSuffix = make_plotTitleSuffix(stimNames);
% groupShading = stimGroupShading;

    % GROUP CHRONOLOGICALLY BY BLOCKS
    trialGroups = zeros(1, nTrials);
    for iBound = 1:numel(groupBounds)-1
        trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
    end
    trialGroups(groupBounds(end)+1:end) = iBound + 1;
    trialGroups = trialGroups .* goodTrials;
    plotTitleSuffix = '';
    fileNameSuffix = '_Blocks_Separated';
    if ~isempty(blockShading)
        sepBlockStims = 1;
    end
    binNames = blockNames;
    plotTitleSuffix = make_plotTitleSuffix(binNames);
    groupShading = blockShading;
%     
% trialGroups = trialGroups - 1;
% trialGroups(1:20) = 0;
% binNames = binNames(2:4);
% 
% % SPLIT TRIALS INTO HIGH-MOVEMENT AND LOW-MOVEMENT GROUPS
% sdCap = 5;
% analysisWindow = [1 1];
% groupFraction = 8;
% offsetAlign = 0;
% 
% % Extract FicTrac speed data
% rawData = infoStruct.ftData.moveSpeed;         % --> [frame, trial]
% rawData = rawData';                            % --> [trial, frame]
% plotData = (rad2deg(rawData .* FRAME_RATE));   % --> [trial, frame] (deg/sec)
% 
% % Cap values at n SD above mean
% capVal = mean(plotData(:), 'omitnan') + (sdCap * std(plotData(:), 'omitnan'));
% plotData(plotData > capVal) = capVal;
% 
% % Smooth data
% smSpeedData = smoothdata(plotData, 2, 'gaussian', smWin);
% 
% stimTrialGroups = stimTrialGroups .* goodTrials;
% groupNums = 1:numel(unique(stimTrialGroups(stimTrialGroups > 0))) * 2;
% trialGroups = []; newGroupShading = []; B = [];
% for iGroup = 1:numel(unique(stimTrialGroups(stimTrialGroups > 0)))
%    
%     currGroupTrials = stimTrialGroups == iGroup;
%     currGroupShading = stimGroupShading{iGroup};
%     currStimStartTime = stimEpochs(currGroupShading, 1 + offsetAlign);
%     
%     % Get mean ficTrac moveSpeed around each stim onset
%     startFrame = round((currStimStartTime - analysisWindow(1)) * FRAME_RATE);
%     endFrame = round((currStimStartTime + analysisWindow(2)) * FRAME_RATE);
%     analysisMean = mean(smSpeedData(currGroupTrials, startFrame:endFrame), 2);
%     
%     [B{iGroup}, I] = sort(analysisMean);
%     
%     nGroupTrials = round(numel(I)/groupFraction);
%     lowMoveTrials = I(1:nGroupTrials);
%     highMoveTrials = I(end - nGroupTrials + 1:end);
%     
%     newGroup = zeros(numel(I), 1);
%     newGroup(lowMoveTrials) = groupNums(1); groupNums(1) = [];
%     newGroup(highMoveTrials) = groupNums(1); groupNums(1) = [];
%     
%     trialGroups(currGroupTrials) = newGroup;
%     newGroupShading = [newGroupShading, {currGroupShading}, {currGroupShading}];
% end
% groupShading = newGroupShading;
% fileNameSuffix = '_HiLowMove_test_plots';
% % plotTitleSuffix = make_plotTitleSuffix({});
% figure(100);clf;hold on; 
% plot(B{1}); plot(B{2}); plot(1:numel(B{1}), zeros(1, numel(B{1})), '--')
% yL = ylim();
% plot([nGroupTrials, nGroupTrials], yL);
% plot([numel(B{1}) - nGroupTrials + 1, numel(B{1}) - nGroupTrials + 1], yL);
% title('Mean speed around stim onset');
% ylabel('Mean speed (mm/sec)')
% legend 1 2

% % 
% % SPLIT TRIALS INTO GROUPS BASED ON BEHAVIORAL RESPONSE TO DIFFERENT STIMULI
% sdCap = 5;
% analysisWindow = [1 2];
% groupFraction = 8;
% offsetAlign = 0;
% 
% % Extract FicTrac speed data
% rawData = infoStruct.ftData.moveSpeed;         % --> [frame, trial]
% rawData = rawData';                            % --> [trial, frame]
% plotData = (rad2deg(rawData .* FRAME_RATE));   % --> [trial, frame] (deg/sec)
% 
% % Cap values at n SD above mean
% capVal = mean(plotData(:), 'omitnan') + (sdCap * std(plotData(:), 'omitnan'));
% plotData(plotData > capVal) = capVal;
% 
% % Smooth data
% smSpeedData = smoothdata(plotData, 2, 'gaussian', smWin);
% 
% stimTrialGroups = stimTrialGroups .* goodTrials;
% groupNums = 1:numel(unique(stimTrialGroups(stimTrialGroups > 0))) * 3;
% trialGroups = []; newGroupShading = []; B = [];
% for iGroup = 1:numel(unique(stimTrialGroups(stimTrialGroups > 0)))
%    
%     currGroupTrials = stimTrialGroups == iGroup;
%     currGroupShading = stimGroupShading{iGroup};
%     currStimStartTime = stimEpochs(currGroupShading, 1 + offsetAlign);
%     
%     % Get mean ficTrac moveSpeed from before and after each stim onset
%     startFrame = round((currStimStartTime - analysisWindow(1)) * FRAME_RATE);
%     endFrame = round((currStimStartTime + analysisWindow(2)) * FRAME_RATE);
%     stimFrame = round(currStimStartTime * FRAME_RATE);
%     baselineMean = mean(smSpeedData(currGroupTrials, startFrame:stimFrame), 2);
%     respMean = mean(smSpeedData(currGroupTrials, stimFrame:endFrame), 2);
% %     respMean = mean(smSpeedData(currGroupTrials, stimFrame+FRAME_RATE:endFrame), 2);
%     
%     meanDiff = (respMean - baselineMean);
%     
%     [B{iGroup}, I] = sort(meanDiff);
%     
%     nGroupTrials = round(numel(I)/groupFraction);
%     lowMoveTrials = I(1:nGroupTrials);
%     highMoveTrials = I(end - nGroupTrials + 1:end);
% %     lowDiffTrials = I(round(numel(I)/2) - round(nGroupTrials/2):round(numel(I)/2) + round(nGroupTrials/2))
%     
%     newGroup = zeros(numel(I), 1);
%     newGroup(lowMoveTrials) = groupNums(1); groupNums(1) = [];
%     newGroup(highMoveTrials) = groupNums(1); groupNums(1) = [];
% %     newGroup(lowDiffTrials) = groupNums(1); groupNums(1) = [];
%     trialGroups(currGroupTrials) = newGroup;
%     newGroupShading = [newGroupShading, {currGroupShading}, {currGroupShading}];
% %     newGroupShading = [newGroupShading, {currGroupShading}, {currGroupShading}, {currGroupShading}];
% end
% groupShading = newGroupShading;
% 
% figure(100);clf;hold on; 
% plot(B{1}); plot(B{2}); plot(1:numel(B{1}), zeros(1, numel(B{1})), '--')
% yL = ylim();
% plot([nGroupTrials, nGroupTrials], yL);
% plot([numel(B{1}) - nGroupTrials + 1, numel(B{1}) - nGroupTrials + 1], yL);
% title('Mean speed increase at stim onset');
% ylabel('Delta speed (mm/sec)')
% legend 1 2

% --------------------------------------------------------------------------------------------------

% % GROUP BY BLOCK TYPE
% groupNames = blockNames(1:2);
% trialGroups = zeros(1, nTrials);
% trialGroups = trialGroups .* goodTrials;
% for iBound = 1:numel(groupBounds)-1
%     trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(logical(mod(trialGroups, 2))) = 1;
% trialGroups(~logical(mod(trialGroups, 2))) = 2;
% fileNameSuffix = '_Block_Type_Groups';
% groupShading 

% --------------------------------------------------------------------------------------------------

try
% Create base file name and select data
if useDff
    fileNamePrefix = ['Behavior_', ftVar, '_Coded_Trial_Avg_Dff_', ...
            regexprep(expDate, {'_', 'exp'}, {'', '-'})];
    flData = ROIDffAvg;
    baselineData = ROIDffAvg(:,:,end);
else
    fileNamePrefix = ['Behavior_', ftVar, '_Coded_Trial_Avg_RawF_', ...
            regexprep(expDate, {'_', 'exp'}, {'', '-'})];
    flData = ROIDataAvg;
    baselineData = ROIDataAvg(:,:,end);
end

% Extract correct FicTrac variable
switch ftVar
    case 'moveSpeed'
        behavData = ftData.(ftVar) * 4.5 * FRAME_RATE;
        cbLabel = 'FicTrac move speed (mm/sec)';
    case 'yawVel'
        behavData = rad2deg(ftData.yawSpeed * pi * FRAME_RATE);
        cbLabel = 'FicTrac yaw velocity (deg/sec)';
    case 'yawSpeed'
        behavData = abs(rad2deg(ftData.yawSpeed * pi * FRAME_RATE));
        cbLabel = 'FicTrac unsigned yaw speed (deg/sec)';
end

% Subtract baseline if necessary
if subtractBaseline
    currNumROIs = nROIs - 1;
    fileNamePrefix = [fileNamePrefix, '_BaselineSub'];
    flData = flData - repmat(baselineData, 1, 1, size(flData, 3));
    flData(:,:,end) = [];
else
    currNumROIs = nROIs;
end

ROIlist = 1:currNumROIs;
%     ROIlist = [1 2 3];
yL = []; f = [];
for iROI = ROIlist       

    nGroups = numel(unique(trialGroups(trialGroups > 0)));
    nRows = nGroups + 1;

    % Update plot title
    if useDff
        plotTitlePrefix = [regexprep(expDate, {'_', 'exp'}, {'', '-'}), ...
                '    ', ROInames{iROI}, '   ', ftVar, '-coded'];
        yAxisLabel = 'dF/F';
    else
        plotTitlePrefix = [regexprep(expDate, {'_', 'exp'}, {'', '-'}), ...
                 '    ', ROInames{iROI}, '   ', ftVar, '-coded'];
        yAxisLabel = 'Raw F';
    end
    if subtractBaseline
        plotTitlePrefix = [plotTitlePrefix, '  (baseline-sub)'];
    end

    % Create figure
    f = figure(iROI); clf
    f.Position = figPos;
    f.Color = [1 1 1];
    
    % Plot reference image and ROI outline
    currPlane = ROImetadata{iROI}(1).plane;
    xData = ROImetadata{iROI}(1).xi;
    yData = ROImetadata{iROI}(1).yi;
    subaxis(nRows, 2, [1 2], 'ML', 0.06, 'MR', 0.015, 'MT', 0.05, 'MB', 0.08, 'SH', 0)
    hold on
    imshow(infoStruct.refImg{currPlane}, [0 infoStruct.MAX_INTENSITY]);
    plot(xData, yData, 'Color', 'g');
    title(ROInames{iROI})
    
    clear ax
    for iGroup = 1:nGroups
        currFlData = flData(:, trialGroups == iGroup, iROI); % --> [volume, trial]

        currBehavData = behavData(:, trialGroups == iGroup);
        subaxis(nRows, 2, 1, iGroup+1, 2,1)  %( (iPlot*3 + 1):(iPlot*3 + 3) ))                                                     
        ax(iGroup) = gca;
        ax(iGroup) = plot_behavior_coded_ROI_data(ax(iGroup), currFlData, currBehavData, ...
                   'SingleTrials', singleTrials, ...
                   'SingleTrialAlpha', singleTrialAlpha, ...
                   'StdDevShading', stdDevShading, ...
                   'OutlierSD', outlierSD, ...
                   'VolumeRate', volumeRate, ...
                   'EdgeColorMode', edgeColorMode, ...
                   'Colormap', cm, ...
                   'SmoothWinSize', smWin);
        ax(iGroup).FontSize = fontSize;
        c = colorbar;
        c.Label.String = cbLabel;
        c.FontSize = fontSize;
        ax(iGroup).Position(3) = 0.83; % So the colorbar label fits inside the figure
        cL = ax(iGroup).CLim;
        ax(iGroup).CLim = [cL(1), cL(2) * CLimCap];
        yL{iGroup} = ylim(ax(iGroup));
        if ~useDff
            ylabel('Raw F');
        end
        if iGroup ~= nGroups 
            xlabel([]);
        end
        hold on
        
        % Add stimulus epoch shading
        if ~isempty(groupShading)
            currGroupShading = groupShading{iGroup};
            currStimEpochs = stimEpochs(currGroupShading, :);
            currStimShadingColors = stimShadingColors(currGroupShading);
        else
            currStimEpochs = stimEpochs;
            currStimShadingColors = stimShadingColors;
        end
        for iEpoch = 1:size(currStimEpochs, 1)
            plot_stim_shading(currStimEpochs(iEpoch, :), 'Axes', ax(iGroup), 'Color', ...
                rgb(currStimShadingColors{iEpoch}));
        end
    end
    ax(1).Title.String = [plotTitlePrefix, plotTitleSuffix];

%     % Scale y-axes to match
%     yMin = min([yL{:}]);
%     yMax = max([yL{:}]);
%     for iGroup = 1:nGroups
%         ylim(ax(iGroup), [yMin yMax]);
%     end


    % Save figure ------------------------------------------------------------------------------
    if saveDir
        fileName = [fileNamePrefix, fileNameSuffix, '_', strrep(ROInames{iROI}, ' ', '_')];
        save_figure(f, saveDir, fileName);
    end
    end
catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end

    %% PLOT CORRELATIONS BETWEEN ROIS
try
stimVols = round(volumeRate * stimEpochs);
saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', ...
    num2str(sid), '\Analysis'], 'Select a save directory');

disp(ROInames)

plotTrials = 1:nTrials;
plotROIs = [1:nROIs-2];%[1 3 5];% 
useDff = 1; 
subtractBaseline = 0;
nBins = 50;

statsVertPos = 0.75;
labelFontSize = 11;
statsFontSize = 9;
titleFontSize = 16;

% All frames
plotVols = 1:nVolumes;
% plotTrials = [1:80, 82:nTrials];
fileNameSuffix = '_allVols';     
plotTitleSuffix = '';  
% 
% % Stim period only
% plotVols = stimVols;
% % plotTrials = [1:80, 82:nTrials];
% fileNameSuffix = '_stimVols';
% plotTitleSuffix = ' (stim vols only)';
% nBins = 50;
% 
% % Non-stim period only
% plotVols = [1:stimVols(1), stimVols(2:end)];
% % plotTrials = [1:80, 82:nTrials];
% fileNameSuffix = '_nonStimVols';
% plotTitleSuffix = ' (non-stim vols only)';


try

% Extract required ROI data
currROIs = numel(plotROIs);
if useDff
    ROIData = ROIDffAvg(:,:,plotROIs);
    fileNameSuffix = ['_dFF', fileNameSuffix];
    plotTitleSuffix = ['  —  dF/F', plotTitleSuffix];
    baselineData = repmat(ROIDffAvg(:, :, end), 1, 1, currROIs);
else
    ROIData = ROIDataAvg(:,:,plotROIs);
    fileNameSuffix = ['_RawF', fileNameSuffix];
    plotTitleSuffix = ['  —  raw F', plotTitleSuffix];
    baselineData = repmat(ROIDataAvg(:, :, end), 1, 1, currROIs);
end
if subtractBaseline
    ROIData = ROIData - baselineData;
end
ROIData = ROIData(plotVols, plotTrials, :);

ROIDataRs = reshape(ROIData, [], size(ROIData, 3));

% for i=1:nROIs
%    ROIDataTrim(ROIDataRs(:,i) > ROIub(i), i) = nan; 
% end

[R P RL RU] = corrcoef(ROIDataRs, 'Rows', 'pairwise');
disp(R)
disp(P)
disp((RU - RL))

f = figure(1);clf;
f.Color = [1 1 1];
for iRow = 1:currROIs
    for jCol = 1:currROIs
        
        if iRow < currROIs && jCol > 1 && jCol > iRow 
            subaxis(currROIs-1, currROIs-1, jCol-1, iRow);hold on
            
            histogram2(ROIDataRs(:,iRow), ROIDataRs(:,jCol), nBins, 'DisplayStyle', 'tile')
            
            xlabel(ROInames{iRow});
            ylabel(ROInames{jCol});
            ax = gca();
            ax.XLabel.FontSize = labelFontSize;
            ax.YLabel.FontSize = labelFontSize;
            ax.XTickLabel = [];
            ax.YTickLabel = [];
            ax.Title.FontWeight = 'normal';
            ax.Title.Units = 'normalized';
            ax.Title.HorizontalAlignment = 'left';
            ax.Title.Position = [0.05 statsVertPos 0];
            if P(iRow,jCol) < 0.01
                pValStr = 'p < 0.01';
            else
                pValStr = ['p = ' num2str(round(P(iRow,jCol), 2), 2)];
            end
            ax.Title.String = {['r = ' num2str(R(iRow,jCol), 2)], pValStr};
            ax.Title.FontSize = statsFontSize;
            axis square
            
            if strcmp(regexp(ROInames{iRow}, '(?<=-).*', 'match'), ...
                        regexp(ROInames{jCol}, '(?<=-).*', 'match'))
                ax.XColor = 'red';
                ax.YColor = 'red';
                ax.LineWidth = 1;
            end
        end
        
        hout = suptitle([regexprep(expDate, '\_', '\\_'), '   sid\_', num2str(sid), ...
                '   Frame-by-frame ROI correlations', plotTitleSuffix]);
        hout.FontSize = titleFontSize;
    end
end

% Save figure ------------------------------------------------------------------------------
if saveDir
    fileName = ['ROI_correlation_histograms', fileNameSuffix];
    for iROI = plotROIs
       fileName = [fileName, '_', strrep(ROInames{iROI}, ' ', '_')]; 
    end
    save_figure(f, saveDir, fileName);
end
    
catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end
    %% PLOT CORRELATIONS BETWEEN SPECIFIC ROI PAIRS ONLY
try
stimVols = round(volumeRate * stimEpochs);
saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', ...
    num2str(sid), '\Analysis'], 'Select a save directory');

plotTrials = 1:nTrials;
disp(ROInames)
ROIpairs_X = [1 2];
ROIpairs_Y = [2 3];

useDff = 1; 
subtractBaseline = 0;
nBins = 45;

figPos = [-1900 380 1900 600];
statsVertPos = 0.85;
labelFontSize = 18;
statsFontSize = 14;
titleFontSize = 18;

% All frames
plotVols = 1:nVolumes; 
plotTrials = 1:30;
fileNameSuffix = '';     
plotTitleSuffix = '  (Trials 1-30)';  

% % Stim period only
% plotVols = stimVols;
% fileNameSuffix = '_stimVols';
% plotTitleSuffix = ' (stim vols only)';
% nBins = 60;

% % Non-stim period only
% plotVols = [1:stimVols(1), stimVols(2:end)];
% fileNameSuffix = '_nonStimVols';
% plotTitleSuffix = ' (non-stim vols only)';

try

% Extract required ROI data
if useDff
    ROIData = ROIDffAvg;
    fileNameSuffix = ['_dFF', fileNameSuffix];
    plotTitleSuffix = ['  —  dF/F', plotTitleSuffix];
    baselineData = repmat(ROIDffAvg(:, :, end), 1, 1, nROIs);
else
    ROIData = ROIDataAvg;
    fileNameSuffix = ['_RawF', fileNameSuffix];
    plotTitleSuffix = ['  —  raw F', plotTitleSuffix];
    baselineData = repmat(ROIDataAvg(:, :, end), 1, 1, nROIs);
end
if subtractBaseline
    ROIData = ROIData - baselineData;
end
ROIData = ROIData(plotVols, plotTrials, :);

% Calculate pairwise correlations
ROIpairData = []; P = []; R = []; ROIgroupNames = [];
for iPair = 1:numel(ROIpairs_X)
   currPairData = cat(3, ROIData(:,:,ROIpairs_X(iPair)), ROIData(:, :, ROIpairs_Y(iPair)));
   currDataRs = reshape(currPairData, [], 2);
   [currR, currP, ~, ~] = corrcoef(currDataRs, 'Rows', 'pairwise');
   P(end + 1) = currP(2, 1);
   R(end + 1) = currR(2, 1);
   ROIpairData{end + 1} = currDataRs;
end

% Create figure
f = figure(1);clf;
f.Color = [1 1 1];
f.Position = figPos;
nPairs = numel(ROIpairData);
for iPair = 1:nPairs
    subaxis(1, nPairs, iPair);hold on; 
    histogram2(ROIpairData{iPair}(:,1), ROIpairData{iPair}(:,2), nBins, 'DisplayStyle', 'tile')
    
    xlabel(ROInames{ROIpairs_X(iPair)})
    ylabel(ROInames{ROIpairs_Y(iPair)})
    ax = gca();
    ax.XLabel.FontSize = labelFontSize;
    ax.YLabel.FontSize = labelFontSize;
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    ax.Title.FontWeight = 'normal';
    ax.Title.Units = 'normalized';
    ax.Title.HorizontalAlignment = 'left';
    ax.Title.Position = [0.05 statsVertPos 0];
    if P(iPair) < 0.01
        pValStr = 'p < 0.01';
    else
        pValStr = ['p = ' num2str(round(P(iPair), 2), 2)];
    end
    ax.Title.String = {['r = ' num2str(R(iPair), 2)], pValStr};
    ax.Title.FontSize = statsFontSize;
    axis square
    hout = suptitle([regexprep(expDate, '\_', '\\_'), '   sid\_', num2str(sid), ...
            '   Frame-by-frame ROI correlations', plotTitleSuffix]);
    hout.FontSize = titleFontSize;
end

% Save figure ------------------------------------------------------------------------------
if saveDir
    fileName = ['ROI_paired_correlation_histograms', fileNameSuffix];
    save_figure(f, saveDir, fileName);
end
    
catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end

%% %================================================================================================
%%% FICTRAC DATA PLOTTING
%% %================================================================================================
% s = stimSepTrials;

% tempTrials = goodTrials;
% tempTrials(31:end) = 0;

currStim = goodTrials; %tempTrials% logical(stimTrialGroups);
smWin = 4;
nBins = 50;
thresh = 0.2;
nROIs = size(ROIDataAvg, 3)-1;
currTrial = 2;

try
    
% Pull out good trials from data 
mmSpeedData = ftData.moveSpeed(:,logical(goodTrials .* currStim)) * FRAME_RATE * 4.5;   % --> [frame, trial] (mm/sec)
dHD = rad2deg(ftData.yawSpeed(:,logical(goodTrials .* currStim)) * FRAME_RATE);         % --> [frame, trial] (deg/sec)
dHDSmooth = repeat_smooth(dHD, 2, 'smWin', smWin);                                 
goodFl = ROIDffAvg(:,logical(goodTrials .* currStim),:);
% goodFl = ROIDataAvg(:,logical(goodTrials .* currStim),:);   % --> [volume, trial, ROI]

% Subtract baseline fluorescence
goodFl = goodFl(:,:,1:nROIs) - repmat(goodFl(:,:,end), 1, 1, nROIs);

% Normalize raw fluorescence for each ROI
for iROI = 1:size(goodFl, 3)
    currFl = goodFl(:,:,iROI);
    goodFl(:,:,iROI) = currFl ./ max(currFl(:));
end

goodTrialNums = 1:nTrials;
goodTrialNums(~logical(goodTrials .* currStim)) = [];
goodTrialTypes = infoStruct.trialType;
goodTrialTypes(~logical(goodTrials .* currStim)) = [];

% Scale fluorescence data so min is zero
goodFl = goodFl - min(goodFl(:));

% Calculate correlations between the fly's speed, yaw, and dF/F
allSpeed = as_vector(repeat_smooth(mmSpeedData(volFrames,:), 2, 'smWin', smWin));
allYaw =  as_vector(dHDSmooth(volFrames,:)); 

corrMat = [allSpeed, allYaw];
for iROI = 1:nROIs
    corrMat(:, end + 1) = as_vector(smoothdata(goodFl(:,:,iROI), 1, 'gaussian', smWin)); % --> columns: [speed, yaw, ROI1, ROI2, ...]
end

[R P RU RL] = corrcoef(corrMat);

% Calculate corrCoeffs for each ROI
allTrialCorrCoeffs = [];
for iROI = 1:nROIs
    currCorrMat = corrcoef(corrMat(:, 1), corrMat(:, iROI + 2));
    allTrialCorrCoeffs(iROI) = currCorrMat(1,2);
end
disp(allTrialCorrCoeffs);


% % Calculate corrCoeffs when both speed and fl data have been normalized for each trial
% speedDataSmooth = repeat_smooth(mmSpeedData(volFrames,:), smWin, 1, 1);
% speedDataNorm = speedDataSmooth - min(speedDataSmooth(:)); goodFlNorm = goodFlNorm - min(goodFlNorm(:));
% for iTrial = 1:size(goodFl, 2) 
%     speedDataNorm(:, iTrial) = speedDataSmooth(:, iTrial) ./ max(speedDataSmooth(:, iTrial));
%     for iROI = 1:size(goodFl, 3) 
%         goodFlNorm(:, iTrial, iROI) = (goodFl(:, iTrial, iROI) ./ max(goodFl(:, iTrial, iROI));
%     end
% end
% corrMat = as_vector(speedDataNorm);
% for iROI = 1:nROIs
%     corrMat(:, end + 1) = as_vector(movmean(goodFlNorm(:,:,iROI), smWin, 1)); % --> columns: [speed, ROI1, ROI2, ...]
% end
% normCorrCoeffs = [];
% for iROI = 1:nROIs
%     currCorrMat = corrcoef(corrMat(:, 1), corrMat(:, iROI + 1));
%     normCorrCoeffs(iROI) = currCorrMat(1,2);
% end
% disp(normCorrCoeffs);
% figure(19); clf; plot(corrMat);xlim([0 220]);


% Find best fitting # of lag frames for each ROI (using MoveSpeed);
xCorrs = []; lags = []; lagInds = []; bestLags = [];
for iROI = 1:nROIs
    [xCorrs(:, iROI), lags(:, iROI)] = xcorr(allSpeed, corrMat(:, iROI + 2), 'coeff');
    [~, lagInds(iROI)] = max(xCorrs(:, iROI));
    bestLags(iROI) = lags(lagInds(iROI), iROI);
end

% Calculate lagged correlation coefficients
lagCorrCoeffs = []; lagPs = [];
lagShiftedCorrMat = zeros(numel(allSpeed), nROIs + 1);
lagShiftedCorrMat(:, 1) = allSpeed;
for iROI = 1:nROIs
    currLag = bestLags(iROI);
    if currLag > 0
        moveSpeedLagged = allSpeed(1:end - currLag);
        corrDataLagged = corrMat(1 + currLag:end, iROI + 2);
        lagShiftedCorrMat(1:end - currLag, iROI + 1) = corrMat(1 + currLag:end, iROI + 2);
    else
        moveSpeedLagged = allSpeed(1:end + currLag);
        corrDataLagged = corrMat(1 - currLag:end, iROI + 2);
        lagShiftedCorrMat(1 - currLag:end, iROI + 1) = corrMat(1:end + currLag, iROI + 2);
    end
   disp(size(moveSpeedLagged))
   [R, P, ~, ~] = corrcoef([moveSpeedLagged, corrDataLagged]);
   lagCorrCoeffs(iROI) = R(1, 2);
   lagPs(iROI) = P(1, 2);
end

% Crop ends off of lag-shifted corrMat
lagShiftedCorrMat = lagShiftedCorrMat(max(abs(bestLags)):end - max(abs(bestLags)), :);


% Create lag-shifted matrix for speed-fl plots




% --------------------------------------------------------------------------------------------------
%  Full dataset plots
% --------------------------------------------------------------------------------------------------

plotNum = 1;
if ~exist('figHandles')
    figHandles = [];
end

% Remove frames if move speed falls below threshold
dffDataVec = []; flDataThresh = [];
mmSpeedDataSmooth = repeat_smooth(mmSpeedData(volFrames, :), 2, 'smWin', smWin);
flDataSmooth = repeat_smooth(goodFl, smWin, 1, 2, 'smWin', smWin);
speedThresh = as_vector(mmSpeedDataSmooth > thresh);
mmSpeedDataVec = mmSpeedDataSmooth(:);
speedDataThresh = mmSpeedDataVec(speedThresh);
flDataVec = [];
for iROI = 1:nROIs
    flDataVec(:,iROI) = as_vector(flDataSmooth(:,:,iROI));
    flDataThresh(:,iROI) = flDataVec(speedThresh, iROI);
end

% Plot 2d histograms for dF/F vs move speed
for iROI = 1:nROIs
    if numel(figHandles) >= plotNum && ishandle(figHandles(plotNum))
        figure(figHandles(plotNum));clf;hold on
        plotNum = plotNum + 1;
    else
        f = figure(plotNum);clf;hold on;
        figHandles(plotNum) = f; plotNum = plotNum + 1;
        f.Color = [1 1 1];
        f.Position = [1930 545 500 450];
    end
    ax = gca();
    ax.FontSize = 20;
    histogram2(speedDataThresh, flDataThresh(:,iROI), nBins, 'DisplayStyle', 'tile')
    title(['ROI ', num2str(iROI), '  (r = ', sprintf('%.3f)', allTrialCorrCoeffs(iROI))]); ...
            xlabel('Move speed (mm/sec)'); ylabel('Raw F (AU)')
    ax.FontSize = 14;
end

% 2D histogram of move speed vs yaw speed
if numel(figHandles) >= plotNum && ishandle(figHandles(plotNum))
    figure(figHandles(plotNum));clf;hold on
    plotNum = plotNum + 1;
else
    f = figure(plotNum);clf;hold on;
    figHandles(plotNum) = f; plotNum = plotNum + 1;
    f.Color = [1 1 1];
    f.Position = [1930 0 500 450];
end
ax = gca();
ax.FontSize = 20;
title('XY vs. Yaw speed')
histogram2(allYaw, smoothdata(as_vector(mmSpeedData(volFrames, :)), 'gaussian', 3), nBins, 'DisplayStyle', 'tile')
xlabel('Yaw speed (deg/sec)'); ylabel('Move speed (mm/sec)')
ax.FontSize = 14;

% Plot ROIs on reference images
if numel(figHandles) >= plotNum && ishandle(figHandles(plotNum))
    figure(figHandles(plotNum));clf;hold on
    plotNum = plotNum + 1;
else
    f = figure(plotNum);clf;hold on;
    figHandles(plotNum) = f; plotNum = plotNum + 1;
    f.Color = [1 1 1];
    f.Position = [2435 0 550 650];
end
for iROI = 1:nROIs
    subaxis(nROIs, 1, iROI, 'M', 0.01)
    ax = gca(); hold on
    ax.FontSize = 20;
    title(['ROI ', num2str(iROI)])
    imshow(refImg{ROImetadata{iROI}(1).plane}, [0 infoStruct.MAX_INTENSITY]);hold on
    plot(ROImetadata{iROI}(1).xi, ROImetadata{iROI}(1).yi, 'LineWidth', 2)
    set(gca, 'FontSize', 14)
end

% Plot raw fluorescence heatmap
for iROI = 1:nROIs
    if numel(figHandles) >= plotNum && ishandle(figHandles(plotNum))
        figure(figHandles(plotNum));clf
        plotNum = plotNum + 1;
    else
        f = figure(plotNum);clf;
        figHandles(plotNum) = f; plotNum = plotNum + 1;
        f.Color = [1 1 1];
        f.Position = [510 545 600 450];
    end
    imagesc(goodFl(:,:,iROI)');
    ax = gca();
    ax.FontSize = 14;
    ax.YTick = 0:10:nTrials;
    ax.XTick = [0:5:trialDuration] * volumeRate;
    ax.XTickLabel = 0:5:trialDuration;
    title(['ROI ', num2str(iROI), '  -  Raw fluorescence (AU)']);
    xlabel('Time (sec)')
    ylabel('Trial')
    
end%iROI

% Plot movement speed heatmap (capping values at 3 SD above mean)
if numel(figHandles) >= plotNum && ishandle(figHandles(plotNum))
    figure(figHandles(plotNum));clf
    plotNum = plotNum + 1;
else
    f = figure(plotNum);clf;
    figHandles(plotNum) = f; plotNum = plotNum + 1;
    f.Color = [1 1 1];
    f.Position = [1115 545 600 450];
end
ax = gca;
capVal = mean(mmSpeedData(:)) + (3 * std(mmSpeedData(:)));
speedDataTrim = mmSpeedData;
speedDataTrim(speedDataTrim > capVal) = capVal;
imagesc(smoothdata(speedDataTrim', 2, 'smWin', smWin));
ax.FontSize = 14;
title('Move speed (mm/sec)');
xlabel('Time (sec)')
ax.XTick = [0:5:trialDuration] * FRAME_RATE;
ax.XTickLabel = 0:5:trialDuration;
ylabel('Trial')
% 
% % Plot closed loop stim
% f = figure(plotNum); clf; plotNum = plotNum + 1;
% f.Color = [1 1 1];
% f.Position = [1115 545 600 450];
% imagesc(annotationTypes{contains(annotationTypeSummary.AnnotationType, 'ClosedLoop')}.frameAnnotArr)
% set(gca, 'FontSize', 14)
% title('Stim delivery');
% xlabel('Volume')
% ylabel('Trial')
% xTickFrame = [0:1/trialDuration:1] * (trialDuration * FRAME_RATE);
% xTickLabels = [0:1/trialDuration:1] * trialDuration;
% ax = gca()
% ax.XTick = xTickFrame;
% ax.XTickLabel = xTickLabels;

% ---------------------------------------------------------------------------------------------------
%  Single trial plots
% ---------------------------------------------------------------------------------------------------
currGroup =  find(currTrial <= groupBounds, 1) - 1;

shadeTimes = stimEpochs(blockShading{(currGroup)},:);
shadeColors = stimShadingColors(blockShading{(currGroup)});
shadeVols = round(shadeTimes * volumeRate);

xTickVol = [0:1/trialDuration:1] * (trialDuration * volumeRate);
xTickLabels = [0:1/trialDuration:1] * trialDuration;

% Plot move speed + fluorescence for current trial
if numel(figHandles) >= plotNum && ishandle(figHandles(plotNum))
    figure(figHandles(plotNum));clf;hold on
    plotNum = plotNum + 1;
    tightVal = 0;
else
    f = figure(plotNum);clf;hold on;
    figHandles(plotNum) = f; plotNum = plotNum + 1;
    f.Color = [1 1 1];
    f.Position = [5 50 950 400];
    f.UserData = f.Position;
end
ax = gca();
ax.FontSize = 20;
ax.XTick = xTickVol;
ax.XTickLabels = xTickLabels;
plotColors = default_plot_colors();
smSpeedData = repeat_smooth(mmSpeedData, 2, 'smWin', smWin);
currMoveSpeed = repeat_smooth(smSpeedData(volFrames,currTrial), 2, 'smWin', smWin);

yyaxis left

plot(currMoveSpeed, 'LineWidth', 2, 'Color', plotColors(1,:))
ylabel('Move speed (mm/sec)');
lgdStr = {'Move speed'};
% ylim([0 max(as_vector(repeat_smooth(mmSpeedData, smWin, 1, 1)))]);
% ylim([0 max(repeat_smooth(currMoveSpeed, smWin, 1, 2))]);
% ylim([0 30])
yyaxis right
for iROI = 1:nROIs
    
    currData = goodFl(:,currTrial, iROI);
    currDataSm = smoothdata(currData, 'gaussian', smWin);
    currDataZeroed = currDataSm - min(currDataSm(:));
    currDataNorm = currDataZeroed ./ max(currDataZeroed(:));
    plot(currDataNorm, ':', 'LineWidth', 2, 'Color', plotColors(iROI + 1, :))
    
%     plot(smooth(goodFl(:,currTrial,iROI),smWin), ':', 'LineWidth', 2, 'Color', plotColors(iROI + 1, :))
    lgdStr{end + 1} = ROInames{iROI};
end
ylabel('Raw florescence (AU)');
legend(lgdStr, 'autoupdate', 'off', 'location', 'nw');
if ~isempty(infoStruct.stimOnsetTimes)
    for iEpoch = 1:size(shadeVols, 1)
        plot_stim_shading(shadeVols(iEpoch, :), 'Axes', ax, 'Color', rgb(shadeColors{iEpoch}));
    end
else
    annotNames = annotationTypeSummary.AnnotationType;
    if any(contains(annotNames, 'ClosedLoopStim'))
        clAnnot = annotationTypes{contains(annotNames, 'ClosedLoopStim')};
        for iEvent = 1:clAnnot.nEvents 
            evList = clAnnot.eventList;
            if evList(iEvent, 3) == currTrial
                plot_stim_shading([evList(iEvent,1), evList(iEvent,2)], 'Axes', ax)                
            end
        end       
    end
end
title(['Move speed  -  Trial #', num2str(goodTrialNums(currTrial)), ' - ', goodTrialTypes{currTrial}])
xlabel('Time (sec)')
% ylim([0 max(goodFl(:))])
% ylim([0 max(as_vector(goodFl(:, currTrial, :)))])
% ylim([0 max(goodFl(:)) * 0.6])
tight_ax(ax);
ax.FontSize = 14;

% Plot normalized yaw speed + dF/F for current trial
if numel(figHandles) >= plotNum && ishandle(figHandles(plotNum))
    figure(figHandles(plotNum));clf;hold on
    plotNum = plotNum + 1;
else
    f = figure(plotNum);clf;hold on;
    figHandles(plotNum) = f; plotNum = plotNum + 1;
    f.Color = [1 1 1];
    f.Position = [960 50 950 400];
    f.UserData = f.Position;
end
ax = gca();
ax.FontSize = 20;
ax.XTick = xTickVol;
ax.XTickLabels = xTickLabels;
currHD = dHDSmooth(volFrames,currTrial);
yyaxis left
plot(repeat_smooth(currHD, 2, 'smWin', smWin), 'LineWidth', 2, 'Color', plotColors(1,:))
ylabel('Yaw speed (deg/sec)');
lgdStr = {'Yaw speed'};
yyaxis right
for iROI = 1:nROIs
    
%      currData = goodFl(:,currTrial, iROI);
%     currDataSm = smooth(currData, smWin);
%     currDataZeroed = currDataSm - min(currDataSm(:));
%     currDataNorm = currDataZeroed ./ max(currDataZeroed(:));
%     plot(currDataNorm, ':', 'LineWidth', 2, 'Color', plotColors(iROI + 1, :))
%     
    
    
    
    plot(smoothdata(goodFl(:,currTrial,iROI), 'gaussian', smWin), '--', 'LineWidth', 2, 'Color', plotColors(iROI + 1, :))
    lgdStr{end + 1} = ROInames{iROI};
end
ylabel('Raw florescence (AU)');
legend(lgdStr, 'autoupdate', 'off', 'location', 'nw');
if ~isempty(infoStruct.stimOnsetTimes)
    for iEpoch = 1:size(shadeVols, 1)
        plot_stim_shading(shadeVols(iEpoch, :), 'Axes', ax, 'Color', rgb(shadeColors{iEpoch}));
    end
else
    annotNames = annotationTypeSummary.AnnotationType;
    if any(contains(annotNames, 'ClosedLoopStim'))
        clAnnot = annotationTypes{contains(annotNames, 'ClosedLoopStim')};
        for iEvent = 1:clAnnot.nEvents 
            evList = clAnnot.eventList;
            if evList(iEvent, 3) == currTrial
                plot_stim_shading([evList(iEvent,1), evList(iEvent,2)], 'Axes', ax)                
            end
        end       
    end
end
title(['Yaw speed  -  Trial #', num2str(goodTrialNums(currTrial)), ' - ', goodTrialTypes{currTrial}])
xlabel('Time (sec)')
tight_ax(ax);
ax.FontSize = 14;

% Scatter plots of current trial dF/F vs move speed
for iROI = 1:nROIs
    if numel(figHandles) >= plotNum && ishandle(figHandles(plotNum))
        figure(figHandles(plotNum));clf;hold on
        plotNum = plotNum + 1;
    else
        f = figure(plotNum);clf;hold on;
        figHandles(plotNum) = f; plotNum = plotNum + 1;
        f.Color = [1 1 1];
        f.Position = [5 545 500 450];
    end
    ax = gca();
    ax.FontSize = 14;
    currMoveSpeed = mmSpeedData(volFrames,currTrial);
    plot(smoothdata(currMoveSpeed, 'gaussian', smWin), smoothdata(goodFl(:, currTrial, iROI), 'gaussian', smWin), 'o', 'color', 'b', 'markersize', 5);
    C = corrcoef(smoothdata(currMoveSpeed, 'gaussian', smWin), smoothdata(goodFl(:, currTrial, iROI), 'gaussian', smWin));
    title(['ROI ', num2str(iROI),'  (r = ', sprintf('%.3f)', C(2,1))]);
    tight_ax(ax);
    xlabel('Move speed (mm/sec)'); ylabel('Raw F (AU)')
end

catch foldME; rethrow(foldME); end
%% =================================================================================================
%  Create vid with trial-by-trial FicTrac plots
%% =================================================================================================
try
saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');

% Extract FicTrac variables
mmXY = ftData.intXY * 4.5; % --> [frame, var, trial] (mm)
HD = ftData.intHD;    % --> [frame, trial]
mmSpeedData = ftData.moveSpeed * FRAME_RATE * 4.5; % --> [frame, trial] (mm/sec)
dHD = rad2deg(ftData.yawSpeed * FRAME_RATE);
rate = 2 * (12/25);
[kb, ka] = butter(2,rate); % LP filter for velocity data

% % Create VidWriter
% myVidWriter = VideoWriter(fullfile(saveDir, 'FicTrac_Trial_Plots.avi'));
% myVidWriter.FrameRate = 1;
% open(myVidWriter)

for iTrial = 1:size(mmXY, 3)
    
    disp(num2str(iTrial))
    
    currIntXY = mmXY(:, :, iTrial);
    currFrameCount = ftData.frameCounter(:, iTrial);
    currSeqNum = ftData.seqNum(:, iTrial);
    currMoveSpeed = mmSpeedData(:, iTrial);
    currHD = HD(:, iTrial);
    currYawSpeed = dHD(:, iTrial);
    currForwardMove = ftData.intForwardMove(:, iTrial);
    currSideMove = ftData.intSideMove(:, iTrial);

    % Smooth velocity data
    smoothWin = 5;
    smoothVelData = repeat_smooth(currMoveSpeed, 2, 'smWin', smoothWin);  % --> [sample, axis, trial]

%     % LP filter velocity data
%     velData = filtfilt(kb, ka, smoothVelData);
    
    % Heading direction data
    uwHD = unwrap(currHD);
    
    % Create figure
    h = figure(100); clf
    screenSize = get( groot, 'Screensize' );
    h.Color = [1 1 1];
    h.OuterPosition = [50 100 (screenSize(3) - 100) screenSize(4) - 150];
    
    % Create axes
    M = 0.01;
    P = 0.00;
    axMove = subaxis(3,6,[1:3 7:9], 'S', 0, 'M', M, 'PB', 0.05, 'PT', 0.05);
    axAnnot = subaxis(3,6,[4:6], 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.00);
    axVel = subaxis(3,6,[10:12], 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.00);
    axYawSpeed = subaxis(3,6,[16:18], 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.00);
    axHD = subaxis(3,6,[13:15], 'S', 0, 'M', M, 'PB', 0.05, 'PR', 0.02);
    
    % Movement map plot
    axes(axMove);
    hold on
    nSegs = trialDuration;
    vectorLen = floor(size(currIntXY, 1) / nSegs);
    cm = jet(nSegs);
    for iSeg = 1:nSegs
        if iSeg == 1
            pad = 1;
        else
            pad = 0;
        end
        currSegInd = (pad + ((iSeg-1)*vectorLen)):(iSeg * vectorLen);
        plot(axMove, smoothdata(currIntXY(currSegInd, 1), 'gaussian', 3), smoothdata(currIntXY(currSegInd, 2), 'gaussian', 3), 'Color', cm(iSeg, :), 'LineWidth', 2)
    end
    axis equal;
    lims = 1.1 * (max(abs([ylim(axMove), xlim(axMove)])));
    if lims < 1
        lims = 1;
    end
    xlim([-lims lims])
    ylim([-lims lims])
    legend({'2D movement (mm)'}, 'FontSize', 11)
    title([infoStruct.trialType{iTrial}, ' (', num2str(infoStruct.stimOnsetTimes(iTrial)), ...
          '-', num2str(infoStruct.stimOnsetTimes(iTrial) + infoStruct.stimDurs(iTrial)), ' sec)'])
    x = currIntXY(iFrame, 1);
    y = currIntXY(iFrame, 2);
    
    % Head direction plot
    axes(axHD)
    hold on
    plot(uwHD, 'linewidth', 2, 'Color', 'b');
    legend({'Heading (rad)'}, 'FontSize', 11, 'Location', 'nw', 'AutoUpdate', 'off');
    axHD.XTick = xTickFR;
    axHD.XTickLabel = xTickLabels;
    xlabel('Time (sec)')
    plot_stim_shading(shadeFrames, 'Axes', axHD);
    
    % Yaw speed plot
    axes(axYawSpeed);
    hold on
    plot(smoothdata(currYawSpeed, 'gaussian', 11), 'LineWidth', 2, 'color', 'r');
    axYawSpeed.XTick = xTickFR;
    axYawSpeed.XTickLabel = xTickLabels;
    legend({'Yaw speed (deg/sec)'}, 'FontSize', 11, 'Location', 'nw', 'AutoUpdate', 'off')
    xlabel('Time (sec)')
    plot_stim_shading(shadeFrames, 'Axes', axYawSpeed);
    
    % Velocity plot
    axes(axVel)
    hold on
    plot(smoothVelData, 'LineWidth', 2, 'Color', rgb('green'));
    axVel.XTick = xTickFR;
    axVel.XTickLabel = xTickLabels;
    legend({'Speed (mm/sec)'}, 'FontSize', 11, 'AutoUpdate', 'off', 'location', 'nw')
    plot_stim_shading(shadeFrames, 'Axes', axVel);
    
    % Anvil annotations
    annotArr = behaviorAnnotArr;
    annotArr(annotArr > 0) = annotArr(annotArr > 0) - 1;
    axes(axAnnot)
    hold on
    plot(annotArr(iTrial, :), '-*', 'Color', 'k')
    cm = [rgb('Navy'); rgb('Cyan')*0.9; rgb('green'); rgb('Gold')];
    
    for iPlot = 1:4
        currData = annotArr(iTrial, :);
        currData(currData ~= iPlot - 1) = nan;
        plot(currData, '*', 'color', cm(iPlot, :))
    end
    ylim([-1 4])
    axAnnot.YTick = [0 1 2 3];
    axAnnot.YTickLabel = {'Quiescence', 'Locomotion', 'Grooming', 'IsolatedMovement'};
    axAnnot.XTick = xTickFR;
    axAnnot.XTickLabel = xTickLabels;
    plot_stim_shading(shadeFrames, 'Axes', axAnnot);
     
%     % Write frame to video(s)
%     writeFrame = getframe(h);
%     writeVideo(myVidWriter, writeFrame);
end

% close(myVidWriter)



fs1 = 500;
t1 = 0:1/fs1:1;
fs2 = 129;
t2 = 0:1/fs2:1;

usFl = resample(currFl, 500, 129);

figure(1);clf;hold on
plot(t1(2:end), currFWSpeed/max(currFWSpeed));
plot(t2(2:end), (currFl-min(currFl))/max(currFl-min(currFl)));

catch foldME; rethrow(foldME); end
%% =================================================================================================
%           OTHER ANALYSES                                   
%%%=================================================================================================

    %% PCA
try
% Pull out data for one plane
planeNum = 4;
planeData = pcaData(:,:,:,planeNum);

figure(2); clf;
subplot(2,2,1)
imshow(infoStruct.refImg{planeNum},[0 infoStruct.MAX_INTENSITY])
colormap(gca, 'gray')
% colormap('parula')
for iPlot = 2:4
    subplot(2, 2, iPlot); imagesc(imgaussfilt(planeData(:,:,iPlot-1),0.5));
    colormap(gca, 'bluewhitered')
end

figure(3); clf;
subplot(2,2,1)
for iPlot = 1:4
    subplot(2, 2, iPlot); imagesc(imgaussfilt(planeData(:,:,iPlot+3)));
    colormap(gca, 'bluewhitered')
end

clear planeNum
catch foldME; rethrow(foldME); end
    %% VIEW RAW DATA FOR A SINGLE TRIAL AND PLANE
try
if ~exist('m', 'var')
    [dataFile, sessionDataPath, ~] = uigetfile('*.mat', 'Select an imaging data session file', 'D:\Dropbox (HMS)\2P Data\Imaging Data\');
    if dataFile == 0
        % Skip loading if user clicked "Cancel"
        disp('Initialization cancelled')
        outputMetadata = [];
        dataFileObj = [];
    else
        m = matfile([sessionDataPath, dataFile]); % Only field is 'wholeSession'
    end 
end

planeNum = 9;
trialNum = 1; % Does not account for any skipped trials
    
preview_trial_movie(m, planeNum, trialNum, [], [], []);
clear planeNum trialNum
catch foldME; rethrow(foldME); end
    %% PLOT dF/F WITHIN ROIs THROUGHOUT EXPERIMENT
try
    % Concatate the trial-by-trial dF/F values into one long vector
    nROIs = size(ROIDffAvg, 3);
    concatROIDff = reshape(ROIDffAvg, nVolumes * nTrials, nROIs);   % --> [volume, ROI]
    
    % Get linear arrays of event annotations
    behavVols = annotationTypes{6}.eventVolsLin;
%     odorAVols = annotationTypes{2}.eventVolsLin;
%     odorBVols = annotationTypes{3}.eventVolsLin;
%     noStimVols = annotationTypes{4}.eventVolsLin;
%     soundVols = annotationTypes{5}.eventVolsLin;
%     groomVols = annotationTypes{8}.eventVolsLin;
    odorVols = annotationTypes{1}.eventVolsLin;
    for iROI = 1:nROIs
        
        f = figure(iROI); clf; hold on
        f.Position = [-1850 200 1800 600];
        title(num2str(iROI))
        
        % Plot dF/F and event annotations
        plot(smoothdata(concatROIDff(:,iROI), 'gaussian', 3), 'b');
        plot(odorVols, 'r');
%         plot(odorBVols, 'g');
%         plot(soundVols, 'y');
%         plot(noStimVols, 'm');
        plot(behavVols, 'k');
%         plot(groomVols, 'c');
        
        % Add trial delineators and numbers
        allVols = 1:(nTrials * nVolumes);
        yL = ylim();
        for iVol = allVols(~logical(mod(allVols, nVolumes)))
            plot([iVol, iVol], yL,  'Color', 'k')
            text(iVol + 50, 0.8 * yL(2), num2str(round(iVol / nVolumes)+1));
        end
        
        xlim([0, 3 * nVolumes]); 
        legend({'', 'EtOH','Fly Movements'});
    end
     % 
% % GROUP BY EARLY/LATE FOR A SINGLE STIM TYPE
% targetStim = s.OdorB;
% stimName = 'CO2';
% groupNames = [];
% groupBounds = [1, 40, 80];
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
%    groupNames{iBound} = [stimName, ' trials ', num2str(groupBounds(iBound)), '-', num2str(groupBounds(iBound + 1))];
% end
% trialGroups(groupBounds(end):end) = numel(groupBounds);
% trialGroups(logical(~targetStim .* goodTrials)) = 0;
% groupNames{end + 1} = [stimName, ' trials ', num2str(groupBounds(end)), '-', num2str(nTrials)];
% fileNameSuffix = ['_EarlyVsLate_', stimName];

% % GROUP BY EARLY/LATE FOR TRIALS FOLLOWING TARGET STIM
% targetStim = s.NoStim;
% stimName = 'NoStim';
% groupNames = [];
% groupBounds = [1, 60, 120];
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
%    groupNames{iBound} = ['After ', stimName, ' trials ', num2str(groupBounds(iBound)), '-', num2str(groupBounds(iBound + 1))];
% end
% trialGroups(groupBounds(end):end) = numel(groupBounds);
% trialGroups(logical(~[0, targetStim(1:end-1)] .* goodTrials)) = 0;
% groupNames{end + 1} = ['After ', stimName, ' trials ', num2str(groupBounds(end)), '-', num2str(nTrials)];
% fileNameSuffix = ['_', actionName, '_After_', stimName, '_EarlyVsLate'];
catch foldME; rethrow(foldME); end
    %% CREATE FICTRAC + BEHAVIOR VIDS
try
vidDir = 'D:\Dropbox (HMS)\2P Data\Behavior Vids\2018_06_29\_Movies';
smWin = 5;

vidFiles = dir(fullfile(vidDir, ['sid*tid*.mp4']));

% Pull out good trials from data
mmXYdata = ftData.intXY(:,:,goodTrials) * 4.5;                     % --> [frame, axis, trial] (mm)       
mmSpeedData = ftData.moveSpeed(:,goodTrials) * FRAME_RATE * 4.5;   % --> [frame, trial] (mm/sec)
mmFWSpeed = ftData.fwSpeed(:,goodTrials) * FRAME_RATE * 4.5;       % --> [frame, trial] (mm/sec)
HD = rad2deg(ftData.intHD(:,goodTrials));                          % --> [frame, trial] (deg)
dHD = rad2deg(ftData.yawSpeed(:,goodTrials) * FRAME_RATE);         % --> [frame, trial] (deg/sec)

% goodFl = ROIDffAvg(:,logical(goodTrials .* currStim),:);
goodFl = ROIDataAvg(:,logical(goodTrials),:);          % --> [frame, trial, ROI]
goodTrialNums = 1:nTrials;
goodTrialNums(~goodTrials) = [];
goodTrialTypes = infoStruct.trialType;
goodTrialTypes(~goodTrials) = [];

% Scale fluorescence data so min is zero
goodFl = goodFl - min(goodFl(:));
goodFlNorm = [];
for iROI = 1:nROIs
    goodFlNorm(:,:,iROI) = goodFl(:,:,iROI) / max(as_vector(goodFl(:,:,iROI)));
end

for iTrial = 1:size(goodFl, 2)
    
    vidName = vidFiles(iTrial).name;
    disp(vidName)
    vidFile = fullfile(vidDir, vidName);
    
    % Prepare ficTrac variables
    currX = smoothdata(mmXYdata(:,1,iTrial), 'gaussian', smWin);
    currY = smoothdata(mmXYdata(:,2,iTrial), 'gaussian', smWin);
    currSpeed = smoothdata(mmSpeedData(:,iTrial), 'gaussian', smWin);
    currFWSpeed = smoothdata(mmFWSpeed(:,iTrial), 'gaussian', smWin);
    currYawSpeed = smoothdata(dHD(:,iTrial), 'gaussian', smWin);
    currHD = HD(:,iTrial);
    
    % Calculate shade volumes if applicable
    shadeVols = [];
    if ~isempty(infoStruct.stimOnsetTimes)
        onsetTimes = infoStruct.stimOnsetTimes;
        onsetTimes(~goodTrials) = [];
        stimDurs = infoStruct.stimDurs;
        stimDurs(~goodTrials) = [];
        shadeTimes = [infoStruct.stimOnsetTimes(iTrial), infoStruct.stimOnsetTimes(iTrial) + infoStruct.stimDurs(iTrial)];
        shadeVols = round(shadeTimes * volumeRate);
    end
    
    % Create VidReader and VidWriter
    myVid = VideoReader(vidFile);
    myVidWriter = VideoWriter(fullfile(vidDir, ['Combined_', vidName]), 'MPEG-4');
    myVidWriter.FrameRate = FRAME_RATE;
    open(myVidWriter)
    
    h = figure(10);
    for iFrame = 1:nFrames
        
        currFrame = uint8(readFrame(myVid));
        clf
        ySize = size(currFrame, 1);
        xSize = size(currFrame, 2);
        
        % Create figure
        screenSize = get( groot, 'Screensize' );
        h.Color = [1 1 1];
        h.OuterPosition = [50 100 (screenSize(3) - 150) screenSize(4) - 150];
        xyRatio = h.Position(3) / h.Position(4);
        
        % Create axes
        clf
        M = 0.006;
        P = 0.00;
        axVid = subaxis(2,4, 1, 'S', 0, 'M', M, 'PB', 0.00, 'PR', 0.01);
        axMove = subaxis(2,4, 5, 'S', 0, 'M', M, 'PB', 0.03, 'PR', 0.01);
        axFW = subaxis(2,4, 2:4, 'S', 0, 'M', M, 'PB', 0.06, 'PR', 0.05, 'PL ', 0.03);
        axYaw = subaxis(2,4, 6:8, 'S', 0, 'M', M, 'PB', 0.07, 'PR', 0.05, 'PL ', 0.04);
        
        % Movie frame plot
        imshow(currFrame, [], 'Parent', axVid);
        axVid.XTickLabel = [];
        axVid.Title.String = ['Trial #', num2str(goodTrialNums(iTrial)), '  -  ', goodTrialTypes{iTrial}];
        
        % Movement map plot
        axes(axMove);
        hold on
        nSegs = trialDuration;
        vectorLen = floor(numel(currX) / nSegs);
        cm = jet(nSegs);
        for iSeg = 1:nSegs
            if iSeg == 1
                pad = 1;
            else
                pad = 0;
            end
            currSegInd = (pad + ((iSeg-1)*vectorLen)):(iSeg * vectorLen);
            plot(axMove, currX(currSegInd), currY(currSegInd), 'Color', cm(iSeg, :))
        end
        axis equal; %axis square
        lims = 1.1 * max(abs([ylim(axMove), xlim(axMove)]));
        if lims < 1
            lims = 1;
        end
        xlim([-lims lims])
        ylim([-lims lims])
        legend({'2D movement (mm)'}, 'FontSize', 11, 'Autoupdate', 'off')
        x = currX(iFrame);
        y = currY(iFrame);
        [arrowVec(1), arrowVec(2)] = pol2cart(deg2rad(currHD(iFrame))+(pi * 0), 0.5);
        arrow([x - arrowVec(1)/2, y-arrowVec(2)/2], [x + arrowVec(1)/2, y + arrowVec(2)/2], 'FaceColor', 'green');        
        
        % Plot FW move speed + fluorescence for current trial
        axes(axFW); hold on
        axFW.FontSize = 14;
        plotColors = default_plot_colors();
        fwSpeedSmooth = repeat_smooth(currFWSpeed(:), 2, 'smWin', smWin);
        yyaxis left
        lgdObj = plot(frameTimes, fwSpeedSmooth, 'LineWidth', 2, 'Color', plotColors(1,:));
        plot(frameTimes(iFrame), fwSpeedSmooth(iFrame), 'd', 'markersize', 10, 'Color', 'k', 'markerfacecolor', 'k');
        ylabel('FW move speed (mm/sec)');
        lgdStr = {'FW move speed'};
        ylim([min(currFWSpeed), max(currFWSpeed)]);
        if ~isempty(infoStruct.stimOnsetTimes)
            plot_stim_shading(shadeTimes, 'Axes', axFW);
        end
        yyaxis right
        for iROI = 1:nROIs
            currFl = smoothdata(goodFlNorm(:,iTrial,iROI),'gaussian', smWin);
            lgdObj(end + 1) = plot(volTimes, currFl, ':', 'LineWidth', 2, 'Color', plotColors(iROI + 1, :));
            lgdStr{end + 1} = ['ROI ', num2str(iROI)];
            [~, idx] = min(abs(volFrames - iFrame));
            plot(volTimes(idx), currFl(idx), 'p', 'markersize', 14, 'Color', 'k', 'markerfacecolor', 'k');
        end
        ylabel('Raw florescence (AU)');
        legend(lgdObj, lgdStr, 'autoupdate', 'off');
        ylim([min(squeeze(goodFlNorm(:,iTrial, :))), max(squeeze(goodFlNorm(:, iTrial, :)))])
        axFW.FontSize = 14;
        xlim([0 trialDuration])
        
        
        % Plot yaw speed + fluorescence for current trial
        axes(axYaw); hold on
        axYaw.FontSize = 14;
        plotColors = default_plot_colors();
        yawSpeedSmooth = repeat_smooth(currYawSpeed, 2, 'smWin', smWin);
        yyaxis left
        lgdObj = plot(frameTimes, yawSpeedSmooth, 'LineWidth', 2, 'Color', plotColors(1,:));
        ylabel('Yaw speed (deg/sec)');
        lgdStr = {'Yaw speed'};
        ylim([min(currYawSpeed(:)), max(currYawSpeed(:))]);
        if ~isempty(infoStruct.stimOnsetTimes)
            plot_stim_shading(shadeTimes, 'Axes', axYaw);
        end
        plot(frameTimes(iFrame), yawSpeedSmooth(iFrame), 'd', 'markersize', 10, 'Color', 'k', 'markerfacecolor', 'k');
        yyaxis right
        for iROI = 1:nROIs
            currFl = smoothdata(goodFlNorm(:,iTrial,iROI),'gaussian', smWin);
            lgdObj(end + 1) = plot(volTimes, currFl, ':', 'LineWidth', 2, 'Color', plotColors(iROI + 1, :));
            lgdStr{end + 1} = ['ROI ', num2str(iROI)];
            [~, idx] = min(abs(volFrames - iFrame));
            plot(volTimes(idx), currFl(idx), 'p', 'markersize', 14, 'Color', 'k', 'markerfacecolor', 'k');
        end
        ylabel('Raw florescence (AU)');
        legend(lgdObj, lgdStr, 'autoupdate', 'off');
        xlabel('Time (sec)')
        ylim([min(squeeze(goodFlNorm(:,iTrial, :))), max(squeeze(goodFlNorm(:, iTrial, :)))])
        axYaw.FontSize = 14;
        xlim([0 trialDuration])       
        
        % Write frame to video(s)
        writeFrame = getframe(h);
        writeVideo(myVidWriter, writeFrame);
    end%iFrame
    
    close(myVidWriter)
end
catch foldME; rethrow(foldME); end