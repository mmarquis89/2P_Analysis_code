%% LOAD REGISTERED DATA FILE AND BEHAVIORAL ANNOTATION DATA


[dataFile, parentDir, ~] = uigetfile('*.mat', 'Select an imaging data session file', 'D:\Dropbox (HMS)\2P Data\Imaging Data\');

try
    if dataFile == 0
        % Skip loading if user clicked "Cancel"
        disp('Initialization cancelled')
    else
        % Load matfile object to query session data size
        disp('Loading matfile...')
        m = matfile([parentDir, dataFile]); % Only field is 'wholeSession'
        
        % Load analysis metadata
        disp('Loading analysis metadata...')
        load(fullfile(parentDir, 'analysisMetadata.mat')) % --> 'analysisMetadata' info struct
        infoStruct = analysisMetadata;
        
        % Load reference images file
        disp('Loading reference images...')
%         [refImgFile, refImgPath, ~] = uigetfile('*.mat', 'Select a reference image file', parentDir);
        refImgPath = parentDir; refImgFile = ['refImages_Reg.mat'];
        load(fullfile(refImgPath, refImgFile)) % --> 'refImages', 'channelNum'
        clear refImgPath refImgFile
        
        % Load PCA data
        disp('Loading PCA data...')
        if exist(fullfile(parentDir, ['PCA_data_', dataFile ]), 'file')
            load(fullfile(parentDir, ['PCA_data_', dataFile ])) % --> 'explained', 'pcaData' ([pc, plane], [y, x, pc, plane]
        end
        
        % Load annotation type data
        disp('Loading annotation type data...')
        load(fullfile(parentDir, 'annotationTypes.mat')) % --> 'annotationTypes', 'annotationTypeSummary'
        
        % Load FicTrac data
        if ~isfield(infoStruct, 'ftData')
            disp('Loading FicTrac data...')
            ftDir = fullfile('D:\Dropbox (HMS)\2P Data\Behavior Vids\', infoStruct.expDate, '_Movies\FicTracData');
            if isdir(ftDir)
                ftData = load_fictrac_data(infoStruct, 'Sid', infoStruct.sid, 'ParentDir', ftDir);
            else
                ftData = load_fictrac_data(infoStruct, 'Sid', infoStruct.sid);
            end
            infoStruct.ftData = ftData;
            analysisMetadata = infoStruct;
            save(fullfile(parentDir, 'analysisMetadata.mat'), 'analysisMetadata', '-v7.3');
        else
            disp('FicTrac data already exists in infoStruct...skipping FicTrac loading process')
            ftData = infoStruct.ftData;
        end
        infoStruct.goodTrials(logical(ftData.resets)) = 0; % Omit any trials in which FicTrac reset from analysis
        
        % Load volume-averaged raw session data
        volAvgDataFile = fullfile(parentDir, ['sid_', num2str(infoStruct.sid), '_volAvgSessionData.mat']);
        if exist(volAvgDataFile)
            load(volAvgDataFile) % "volAvgSessionData"
        end
        
        % Reload analysis metadata in case .goodTrials has been updated
        load(fullfile(parentDir, 'analysisMetadata.mat')) % --> 'analysisMetadata' info struct
        infoStruct = analysisMetadata;
        infoStruct.refImg = refImages;
        disp('All data loaded')
        

        
        % ------- Copy variables for convenience -------
        sessionSize = size(m, 'wholeSession');
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
        stimSepTrials = infoStruct.stimSepTrials;
        behaviorAnnotArr = annotationTypes{contains(annotationTypeSummary.AnnotationType, 'move')}.frameAnnotArr;
        
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

clear all; % to make sure variables don't get mixed up with an imaging experiment

sid = 0;
expDir = uigetdir('D:\Dropbox (HMS)\2P Data\Behavior Vids\', 'Select an experiment directory');
expDate = regexp(expDir, '201.*', 'match'); expDate = expDate{:};
savePath = fullfile(expDir, ['sid_', num2str(sid)]);

annotFileName = '_Movies\autoAnnotations.mat';

try
% ----------------  Load stim metadata -------------------------------------------------------------
infoStruct = [];
stimDataFiles = dir(fullfile(expDir, ['metadata*sid_', num2str(sid), '*.mat']));
infoStruct.stimOnsetTimes = []; infoStruct.stimDurs = []; infoStruct.trialType = []; infoStruct.outputData = [];
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
annotData = load(fullfile(expDir, annotFileName)); % variables 'trialAnnotations', 'annotParams', 'ftData', 'flowArr', 'goodTrials', 'behaviorLabels', 'frameInfo'
annotData.nFrames = annotData.frameInfo.nFrames;
annotData.frameTimes = annotData.frameInfo.frameTimes;
annotData.FRAME_RATE = annotData.frameInfo.FRAME_RATE;
infoStruct = setstructfields(infoStruct, annotData);


% ----------------  Load FicTrac data --------------------------------------------------------------
ftData = load_fictrac_data(infoStruct, 'sid', sid, 'ParentDir', fullfile(expDir, '\_Movies\FicTracData'));
infoStruct.ftData = ftData;
infoStruct.goodTrials(logical(ftData.resets)) = 0;

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
%% Remove trials from analysis
try
    
goodTrials = infoStruct.goodTrials;
goodTrials(omitTrials) = 0;

catch foldME; rethrow(foldME); end

%% =================================================================================================
%   BEHAVIOR SUMMARIES                                   
%%%=================================================================================================
 
    %% SET UP PLOTTING PARAMETERS

    stimNames = {'EtOH\_neat', 'Benzaldehyde\_e-1'}; % {'EtOH\_e-2'};%
    stimTrialGroups = [s.OdorA + 2 * s.OdorB]; % [s.OdorA];%
    stimGroupNames = {'OdorA', 'OdorB'}; %{'OdorA'};%
    
    stimEpochs =  [4 6; 8 11];%[6 7; 8 11; 10 11];%[4 5;10 13];%[8 11];% [7 10 ; 10 13];%
    stimShadingColors = {'red', 'green', 'red'}; % {'red', 'green'};%{'green'};%
    stimEpochNames = {'Laser', 'Odor', 'Laser'}; % {'Laser', 'Odor'};%{'Odor'}%
    
    
    % groupBounds = [1:40:nTrials]; groupBounds(2:end) = groupBounds(2:end) - 1;
    % groupBounds = [1, 40:60:nTrials-1]; groupBounds(2:end) = groupBounds(2:end) -1;
    groupBounds = [1, 40];
    groupBounds(1) = 0;
    
    blockNames = {'Baseline', 'Photostim', 'OdorOnly'}%{'OdorOnly', 'BW', 'OdorOnly2', 'FW', 'OdorOnly3', 'BW2'}; %{'Odor Only', 'Photostim'}; %{'OdorOnly', 'FW', 'OdorOnly2', 'BW', 'OdorOnly3', 'FW2'};
    blockShading = {2, [1 2], 2}%{2, [1 2], 2, [2 3], 2, [1 2]}; %{2, [2 3], 2, [1 2], 2, [2 3]};
    
    
    %% PLOT 2-D SUMMARY OF BEHAVIOR DATA ANNOTATIONS
try
saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
sepBlockStims = 0;


trialGroups = [goodTrials];
plotTitleSuffix = '';
fileNameSuffix = '_AllTrials';

% GROUP BY STIM TYPE
% trialGroups = stimTrialGroups; 
% trialGroups = trialGroups .* goodTrials;
% plotTitleSuffix = make_plotTitleSuffix(stimNames); %
% fileNameSuffix = make_fileNameSuffix(stimGroupNames);


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
[~, ax, f] = plot_2D_summary(infoStruct, annotArr, ...
    'plotAxes', ax, ...
    'trialGroups', trialGroups, ...
    'titleStr', titleString, ...
    'sampRate', FRAME_RATE, ...
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
        epochInds = blockShading{iBlock};
        currShadeColors = stimShadingColors(blockShading{iBlock});
        draw_stim_lines(currBlockShading, currShadeColors, ...
                        'plotAxes', ax, 'yLims', blockRanges(iBlock, :));      
    end
    
else
    % Draw one set of lines for the entire plot
    draw_stim_lines(stimEpochs, stimShadingColors, 'plotAxes', ax);
end

% Save figure
if saveDir
    fileName = ['2D_Behavior_Annotation_Summary ', fileNameSuffix, '_', ...
                regexprep(expDate, {'_', 'exp'}, {'', '_'})];
    save_figure(saveDir, fileName);
end%if

catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end
    %% PLOT 1-D SUMMARY OF BEHAVIOR DATA ANNOTATIONS
try
%----- Plot 1D trial-averaged movement data -----
s = stimSepTrials; 
actionNames = {'NA', 'IsoMovement', 'Locomotion'};
saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
actionNum = [3]; % locomotionLabel = 3; noActionLabel = 0; isoMovementLabel = 1;
actionName = actionNames{actionNum};
figTitle = regexprep([expDate, '  —  Fly ', actionName, ' throughout trial'], '_', '\\_');

% ALL TRIALS
trialGroups = [goodTrials];
fileNameSuffix = ['_AllTrials_', actionName];
groupNames = {'All trials'};

% % GROUP BY STIM TYPE
% trialGroups = stimTrialGroups .* goodTrials;
% fileNameSuffix = [make_fileNameSuffix(stimGroupNames), '_', actionName]; 
% groupNames = stimNames;
% % % % 
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

% % PLOT AND COLOR EACH BLOCK SEPARATELY
% groupNames = blockNames;
% trialGroups = zeros(1, nTrials);
% trialGroups = trialGroups .* goodTrials;
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(groupBounds(end)+1:end) = numel(groupBounds);
% fileNameSuffix = '_Blocks_Separated';

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
    cm = parula(nGroups);
    cm = [rgb('blue'); rgb('red'); rgb('green'); rgb('magenta'); rgb('cyan'); rgb('gold'); rgb('lime')];
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
    
    ax.XLim = [20 nFrames-20]; % to improve plot appearance
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
    fileName = regexprep(['1D_Behavior_Annotation_Summary', fileNameSuffix, '_', ...
                        regexprep(expDate, {'_', 'exp'}, {'', '_'})], '_', '\_');
    save_figure(saveDir, fileName);
end%if

catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end
    %% PLOT 2-D SUMMARY OF FICTRAC DATA
try
saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
fontSize = 12;
sepBlockStims = 0;

ftVarName = 'moveSpeed'; % 'moveSpeed', 'fwSpeed', 'yawSpeed', 'yawVel'
sdCap = 2;
smWin = 9;
cmName = @parula;
figTitle = [regexprep(expDate, '_', '\\_'), '  —  FicTrac ', ftVarName];

% % % 
% ALL TRIALS
trialGroups = [];
fileNameSuffix = ['_AllTrials'];
figTitleSuffix = '';

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
subplot(211); hist(as_vector(movmean(plotData, smWin, 2)), 100);
plotData(plotData > capVal) = capVal;
subplot(212);  hist(as_vector(movmean(plotData, smWin, 2)), 100)

% Smooth data
smPlotData = movmean(plotData, smWin, 2);

% Create colormap
cm = cmName(numel(unique(smPlotData)));

% Plot data
titleStr = [figTitle, figTitleSuffix];
fileName = ['2D_FicTrac_', ftVarName '_Summary', fileNameSuffix, '_', ...
            regexprep(expDate, {'_', 'exp'}, {'', '_'})];
[~, ax, f] = plot_2D_summary(infoStruct, smPlotData, ...
                'trialGroups', trialGroups, ...
                'titleStr', titleStr, ...
                'sampRate', FRAME_RATE, ...
                'colormap', cm ...
                );
ax.Title.FontSize = fontSize;

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
        epochInds = blockShading{iBlock};
        currShadeColors = stimShadingColors(blockShading{iBlock});
        draw_stim_lines(currBlockShading, currShadeColors, ...
                        'plotAxes', ax, 'yLims', blockRanges(iBlock, :));      
    end
    
else
    % Draw one set of lines for the entire plot
    draw_stim_lines(stimEpochs, stimShadingColors, 'plotAxes', ax);
end

% Save figure
if saveDir
    fileName = ['2D_FicTrac_', ftVarName, '_Summary', fileNameSuffix, '_', ...
                regexprep(expDate, {'_', 'exp'}, {'', '_'})];
    save_figure(saveDir, fileName);
end

catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end
    %% PLOT 1-D SUMMARY OF FICTRAC DATA
try
    
saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');

includeQuiescence = 1;
if ~includeQuiescence
    fileNameSuffix = 'NoQuiescence_';
else
    fileNameSuffix = '';
end
figTitle = [expDate, '  —  Trial-Averaged FicTrac data'];
trialGroups = goodTrials';
smWin = 11;
cm = [];
% 
% % 
% % ALL TRIALS
% trialGroups = [goodTrials];
% fileNameSuffix = [fileNameSuffix, 'AllTrials'];
% groupNames = {'All trials'};

% % % % % 
% % GROUP BY STIM TYPE
% trialGroups = stimTrialGroups .* goodTrials; 
% fileNameSuffix = [fileNameSuffix, 'StimTypeComparison']; 
% groupNames = stimNames;
% % % 
% % % % 
% % GROUP BY EARLY/LATE
% groupNames = [];
% groupBounds = [1, 30, 60];
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
%    groupNames{iBound} = ['Trials ', num2str(groupBounds(iBound)), '-', num2str(groupBounds(iBound + 1))];
% end
% trialGroups(groupBounds(end):end) = numel(groupBounds);
% trialGroups = trialGroups .* goodTrials;
% groupNames{end + 1} = ['Trials ', num2str(groupBounds(end)), '-', num2str(nTrials)];
% fileNameSuffix = [fileNameSuffix, 'EarlyVsLateTrials'];

% 
% % GROUP BY EARLY/LATE FOR A SINGLE STIM TYPE
% targetStim = s.OdorA;
% stimName = 'OdorA';
% groupNames = [];
% groupBounds = [1, 30, 60];
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
%    groupNames{iBound} = [stimName, ' trials ', num2str(groupBounds(iBound)), '-', num2str(groupBounds(iBound + 1))];
% end
% trialGroups(groupBounds(end):end) = numel(groupBounds);
% trialGroups(logical(~targetStim .* goodTrials)) = 0;
% groupNames{end + 1} = [stimName, ' trials ', num2str(groupBounds(end)), '-', num2str(nTrials)];
% fileNameSuffix = [fileNameSuffix, 'EarlyVsLate_', stimName];

% % PLOT AND COLOR EACH BLOCK SEPARATELY
% groupNames = blockNames;
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(groupBounds(end)+1:end) = numel(groupBounds);
% trialGroups = trialGroups .* goodTrials;
% fileNameSuffix = [fileNameSuffix, 'Blocks_Separated'];

try
    
% Extract relevant data
xTickFR = [0:1/trialDuration:1] * (trialDuration * FRAME_RATE);
xTickLabels = [0:1/trialDuration:1] * trialDuration;
mmSpeedData = ftData.moveSpeed * FRAME_RATE * 4.5;  % --> [frame, trial] (mm/sec)
dHD = abs(rad2deg(ftData.yawSpeed * FRAME_RATE));        % --> [frame, trial] (deg/sec)
fwSpeed = (rad2deg(ftData.yawSpeed * FRAME_RATE));        % --> [frame, trial] (deg/sec)%ftData.fwSpeed * FRAME_RATE * 4.5;        % --> [frame, trial  (mm/sec)
yawVel = rad2deg(ftData.yawSpeed * FRAME_RATE);        % --> [frame, trial] (deg/sec)
nFrames = size(mmSpeedData, 1);

% Create figure
f = figure(3); clf; hold on
f.Position = [-1100 50 900 930];
f.Color = [1 1 1];

% Create axes
M = 0.02;
P = 0.00;
axVel = subaxis(3,1,1, 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.05); hold on
axYawVel = subaxis(3,1,2, 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.06); hold on
axYawSpeed = subaxis(3,1,3, 'S', 0, 'M', M, 'PB', 0.06, 'PL', 0.06); hold on

% Plot data
nGroups = length(unique(trialGroups(trialGroups ~= 0)));
if isempty(cm)
    cm = parula(nGroups);
    cm = [rgb('blue'); rgb('red'); rgb('green'); rgb('magenta'); rgb('cyan'); rgb('gold'); rgb('lime')];
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

    meanSpeed = smooth(mean(currXYSpeed, 2, 'omitnan'), smWin);
    meanYawSpeed = smooth(mean(currYawSpeed, 2, 'omitnan'), smWin);
    meanYawVel = smooth(mean(currYawVel, 2, 'omitnan'), smWin);

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
            plot_stim_shading([stimOnset, stimOffset], 'Color', rgb(stimShadingColors{iStim}), 'Axes', ...
                currAxis);
        end
    end
end

suptitle(regexprep([figTitle, '  —  ', fileNameSuffix], '_', '\\_'));

% Save figure
if saveDir
    fileName = regexprep(['1D_FicTrac_Summary_', fileNameSuffix, '_', ...
                            regexprep(expDate, {'_', 'exp'}, {'', '_'})], '_', '\_');
    save_figure(saveDir, fileName);
end%if

catch foldME; rethrow(foldME); end

catch foldME; rethrow(foldME); end
 %% PLOT OVERLAID 2D MOVEMENT DATA
try
startTime = 9;
plotLen = 10;
limScalar = 0.9;
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
cm = [rgb('blue'); rgb('red'); rgb('green'); rgb('magenta'); rgb('cyan'); rgb('gold'); rgb('lime')];

% % PLOT ALL TRIALS COLORED BY TIME
% surfPlot = 1;
% cm = jet(endFrame - startFrame + 1);
% fileNameSuffix = '_chronological';
% 
% % GROUP BY STIM TYPE
% trialGroups =  [s.OdorA + 2 * s.OdorB] .* goodTrials;
% fileNameSuffix = '_OdorAvsOdorBvsNoStim'; 
% groupNames = stimNames;
% 
% % PLOT AND COLOR EACH BLOCK SEPARATELY
% groupNames = blockNames;
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound)+1:groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(groupBounds(end)+1:end) = numel(groupBounds);
% trialGroups = trialGroups .* goodTrials;
% fileNameSuffix = '_Blocks_Separated';

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
    smoothX = smooth(currXY(:, 1), 7);
    smoothY = smooth(currXY(:, 2), 7);
    
    % Rotate so fly is facing upwards at starting point
    smoothHD = mod(smooth(unwrap(HD(:, iTrial)), 11), (2*pi));
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
    saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
    if saveDir
        % Create filename
        fileName = regexprep(['Movement_Overlay_', num2str(startTime), '-', num2str(startTime + plotLen)...
            '_sec', fileNameSuffix, '_', ...
            regexprep(expDate, {'_', 'exp'}, {'', '_'})], '_', '\_');
        save_figure(saveDir, fileName);
    end
end
catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end

%% =================================================================================================
%   LOAD PROCESSED EVENT DATA 
%% =================================================================================================
try
% Load processed event data
[eventDataFile, eventDataPath, ~] = uigetfile('*.mat', 'Select an event data file', parentDir);
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
    load(fullfile(parentDir, ['CombDffAvg_', eventDataFile])) % --> 'combinedDffAvg', 'allCondSummaries', 'allCondNames', 'combFilterVecs'
else
    disp('No event data file selected - loading cancelled')
end
catch foldME; rethrow(foldME); end
%% =================================================================================================
%   ODOR/SOUND STIM HEATMAPS                                
%%%=================================================================================================
try
    % Show summary again
    odorEventName = 'OdorB';
    eventInd = cellfun(@strcmp, primaryEventNames, repmat({odorEventName}, 1, numel(primaryEventNames)));
    currSummary = allCondSummaries{eventInd};
    currCondNames = allCondNames{eventInd};
    disp(currSummary)
    
    % Calculate absolute max dF/F value across all planes and stim types
    currConds = [2 5];
    makeVid = 1;
    sigma = [0.6];   
    rangeType = 'Max';
    rangeScalar = 0.4;
    saveDir = [];
    fileName = [odorEventName, '_Response_Heatmaps'];

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
    rangeScalar = 0.8;
    makeVid = 1;
    saveDir = [];
    fileName = ['All_frame_dFF_', behaviorNames{actionNum}, '_Heatmaps'];
    titleStr = {['dF/F - ', behaviorNames{actionNum}, ' vs. Quiescence']};
    
    % Load dF/F data
    load(fullfile(parentDir, ['actionDff_', behaviorNames{actionNum}, '.mat'])) % --> 'actionDff'

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
    eventInd = ~cellfun(@isempty, strfind(primaryEventNames, 'locomotion'));
    currSummary = allCondSummaries{eventInd};
    disp(currSummary)
                 
    currCondNames = repmat(allCondNames{eventInd}, 2, 1);
    
    % Calculate absolute max dF/F value across all planes and stim types
    currConds = [3 6];
    sigma = [0.6]; 
    rangeType = 'Max';
    rangeScalar = 0.6;
    makeVid = 1;
    saveDir = [];
    fileName = 'Locomotion_Response_Heatmaps';

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
parentDir = ['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)];
metaDataFileName = 'ROI_metadata.mat';
dffDataFileName = 'ROI_Data_Avg.mat';

% Load metadata
load(fullfile(parentDir, metaDataFileName)); % --> ROImetadata(.mask, .xi, .yi, .plane, .color, .refImg)
infoStruct.ROImetadata = ROImetadata;
infoStruct.nROIs = numel(ROImetadata); nROIs = infoStruct.nROIs;

% Load imaging and dF/F data
load(fullfile(parentDir, dffDataFileName)); % --> 'ROIDataAvg', 'ROIDffAvg', 'ROIDataBaseSub' ([volume, trial, ROI])
disp('ROI data loaded')

% Plot ROIs 
plot_ROIs(ROImetadata);

% Plot mean value in each ROI for each trial
volAvgROIData = squeeze(mean(ROIDataAvg, 1)); % --> [trial, ROI]
figure(nROIs + 1);clf;hold on
legendStr = [];
for iROI = 1:nROIs
    currROIData = volAvgROIData(:, iROI);
    zeroedData = currROIData - min(currROIData(:));
    legendStr{iROI} = num2str(iROI);
    plot(zeroedData);
end
legend(legendStr);

clear metaDataFileName dffDataFileName 
catch foldME; rethrow(foldME); end
    %% PLOT 2D HEATMAPS OF ROI FLUORESCENCE THROUGHOUT EXPERIMENT
try
s = stimSepTrials;
saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');

useDff = 0;
% % 
trialGroups = [];
plotTitleSuffix = '';
fileNameSuffix = '_AllTrials';
% 
% trialGroups = stimTrialGroups; 
% plotTitleSuffix = make_plotTitleSuffix(stimNames); %
% fileNameSuffix = make_fileNameSuffix(stimGroupNames);

%---------------------------------------------------------------------------------------------------

for iROI = 1:nROIs-1
    
    % Create title and save file name
    plotTitle = ['ROI ', num2str(iROI), '  -  Raw fluorescence (AU)', plotTitleSuffix];
    if useDff
        fileNamePrefix = 'ROI_dFF_Summary_';
    else
        fileNamePrefix = 'ROI_Fluorescence_Summary_';
    end
    
    % Select data
    if useDff
        flData = ROIDffAvg;
    else
        flData = ROIDataAvg;
    end
%     if exist('ROIDataBaseSub', 'var')
%         flData = ROIDataBaseSub;
%     end
    
    % Create figure
    f = figure(iROI); clf
    f.Color = [1 1 1];
    f.Position = [-1050 45 1020 950];
    ax = axes();
    
    % Create colormap
    cm = parula(numel(unique(flData(:,goodTrials,iROI))));
    flData(:, ~goodTrials, iROI) = max(as_vector(flData(:,goodTrials,iROI))) + 1;
    if sum(~goodTrials)
        cm = [0 0 0; cm(2:end-1, :); 1 1 1];
    end
    
    % Plot raw fluorescence heatmap
    fileName = [fileNamePrefix, 'ROI_', num2str(iROI), '_', expDate, fileNameSuffix];
    plot_2D_summary(infoStruct, flData(:,:,iROI)', ...
        'plotAxes', ax, ...
        'trialGroups', trialGroups, ...
        'titleStr', [regexprep(expDate, '_', '\\_'), '  -  ROI ', num2str(iROI), '   Raw fluorescence (AU)', plotTitleSuffix], ...
        'saveDir', saveDir, ...
        'colormap', cm, ...
        'fileName', fileName);
    
    % Plot stim times
    hold on
    [nStimEpochs, idx] = max(cellfun(@size, stimShading, repmat({1}, 1, numel(stimShading))));
    for iStim = 1:nStimEpochs
        stimOnsetFrame = stimShading{idx}(iStim, 1) * volumeRate;
        stimOffsetFrame = stimShading{idx}(iStim, 2) * volumeRate;
        plot(ax, [stimOnsetFrame, stimOnsetFrame], ylim(), 'Color', 'k', 'LineWidth', 2)
        plot(ax, [stimOffsetFrame, stimOffsetFrame], ylim(), 'Color', 'k', 'LineWidth', 2)
    end
    
end%iROI
catch foldME; rethrow(foldME); end
    %% PLOT SAME 2D HEATMAPS WITH BASELINE (last ROI) SUBTRACTED
try
s = stimSepTrials;
saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');

useDff = 0;
% % 
% trialGroups = [];
% plotTitleSuffix = '';
% fileNameSuffix = '_AllTrials';
% % 
trialGroups = stimTrialGroups; 
plotTitleSuffix = make_plotTitleSuffix(stimNames); %
fileNameSuffix = make_fileNameSuffix(stimGroupNames);

%---------------------------------------------------------------------------------------------------

for iROI = 1:(nROIs - 1)
    
    % Create title and save file name
    if useDff
        fileNamePrefix = 'ROI_dFF_Summary_BaseSub_';
    else
        fileNamePrefix = 'ROI_Fluorescence_Summary_BaseSub_';
    end
    
    % Select data
    if useDff
        flData = ROIDffAvg(:,:,1:end-1);
        baselineData = ROIDffAvg(:,:,end);
    else
        flData = ROIDataAvg(:,:,1:end-1);
        baselineData = ROIDataAvg(:,:,end);
    end
    
    flData = flData - repmat(baselineData, 1, 1, size(flData, 3));
    
    % Create figure
    f = figure(iROI); clf
    f.Color = [1 1 1];
    f.Position = [-1050 45 1020 950];
    ax = axes();
    
    % Create colormap
    cm = parula(numel(unique(flData(:,goodTrials,iROI))));
    flData(:, ~goodTrials, iROI) = max(as_vector(flData(:,goodTrials,iROI))) + 1;
    if sum(~goodTrials)
        cm = [0 0 0; cm(2:end-1, :); 1 1 1];
    end
    
    % Plot raw fluorescence heatmap
    fileName = [fileNamePrefix, 'ROI_', num2str(iROI), '_', expDate, fileNameSuffix];
    plot_2D_summary(infoStruct, flData(:,:,iROI)', ...
        'plotAxes', ax, ...
        'trialGroups', trialGroups, ...
        'titleStr', [regexprep(expDate, '_', '\\_'), '  -  ROI ', num2str(iROI), '   baseline-sub fluorescence (AU)', plotTitleSuffix], ...
        'saveDir', saveDir, ...
        'colormap', cm, ...
        'fileName', fileName);
    
    % Plot stim times
    hold on
    [nStimEpochs, idx] = max(cellfun(@size, stimShading, repmat({1}, 1, numel(stimShading))));
    for iStim = 1:nStimEpochs
        stimOnsetFrame = stimShading{idx}(iStim, 1) * volumeRate;
        stimOffsetFrame = stimShading{idx}(iStim, 2) * volumeRate;
        plot(ax, [stimOnsetFrame, stimOnsetFrame], ylim(), 'Color', 'k', 'LineWidth', 2)
        plot(ax, [stimOffsetFrame, stimOffsetFrame], ylim(), 'Color', 'k', 'LineWidth', 2)
    end
    
end%iROI
    
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
eventName = 'locomotion';
fileNamePrefix = 'Locomotion_responses_';
shadeDur = 0;
eventInd = contains(primaryEventNames, eventName);
currSummary = allCondSummaries{eventInd};
disp(currSummary)

currConds = [3 6];
currCondNames = allCondNames{eventInd}(currConds);

currDffData = ROIEventDff{eventInd}(currConds); % --> {cond}[volume, event, ROI]

saveDir = 0;
saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');


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
    subaxis(plotPos(1), plotPos(2), 1,'ML', 0.05, 'MR', 0.02, 'MT', 0.05, 'MB', 0.08, 'SV', 0.1, 'SH', 0.05)
    hold on
    imshow(infoStruct.refImg{currPlane}, [0 MAX_INTENSITY]);
    plot(xData, yData, 'Color', 'g');
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
        plot_ROI_data(ax{iPlot}, plotDff, 'EventShading', eventShading, 'VolumeRate', volumeRate, 'VolOffset', volOffset, 'OutlierSD', 5); 
        title(regexprep([expDate, ' - ', currCondNames{iPlot}, ' - ', alignStr{iPlot}], '_', '\\_'))
        yLims(iPlot,:) = ylim(ax{iPlot});
    end
    
    % Scale y-axes to match
    yMin = min(yLims(:));
    yMax = max(yLims(:));
    for iPlot = 1:(nPlots - 1)
       ylim(ax{iPlot}, [yMin yMax]); 
    end
    
    % Save figure -------------------------------------------------------------------------------
    if saveDir
        fileName = [fileNamePrefix, 'ROI_', num2str(iROI), '_', expDate];
        export_fig(fullfile(saveDir, fileName), '-png', f);
        if ~isdir(fullfile(saveDir, 'figFiles'))
            mkdir(fullfile(saveDir, 'figFiles'))
        end
        savefig(f, fullfile(saveDir, 'figFiles', fileName));
    end
end
catch foldME; rethrow(foldME); end
    %% PLOT MEAN ROI dF/F THROUGHOUT TRIAL FOR ONE OR MORE STIM TYPES
try
    saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
    plotTitle = regexprep([expDate, make_plotTitleSuffix(stimNames)], '(?<!\\)_', '\\_');    
    fileNamePrefix = 'Whole_Trial_Responses_';
    s = infoStruct.stimSepTrials;
    filterVecs = [];
    for iStim = 1:numel(stimGroupNames)
        filterVecs = [filterVecs; s.(stimGroupNames{iStim})];
    end
    filterVecs = logical(filterVecs .* repmat(goodTrials, size(filterVecs, 1), 1));
    
    singleTrials = 1;
    singleTrialAlpha = 0.5;
    stdDevShading = 1;
    
    ROIlist = 1:size(ROIDffAvg, 3)-1;
%     ROIlist = [1 2 3];
    
    yL = [];
    for iROI = ROIlist       

        nPlots = size(filterVecs, 1);
        nRows = nPlots + 1;
        
        % Create figure
        f = figure(iROI); clf
        f.Position = [-1450 50 800 900];
        f.Color = [1 1 1];
        
        % Plot reference image and ROI outline
        currPlane = ROImetadata{iROI}(1).plane;
        xData = ROImetadata{iROI}(1).xi;
        yData = ROImetadata{iROI}(1).yi;
        subaxis(nRows, 2, [1 2], 'ML', 0.06, 'MR', 0.015, 'MT', 0.05, 'MB', 0.08, 'SH', 0)
        hold on
        imshow(infoStruct.refImg{currPlane}, [0 infoStruct.MAX_INTENSITY]);
        plot(xData, yData, 'Color', 'g');

        clear ax
        for iPlot = 1:nPlots
            currDffAvg = ROIDffAvg(:, filterVecs(iPlot, :), iROI); % --> [volume, trial]
%             currDffAvg = ROIDataAvg(:, filterVecs(iPlot, :), iROI); % --> [volume, trial]
            subaxis(nRows, 2, 1, iPlot+1, 2,1)  %( (iPlot*3 + 1):(iPlot*3 + 3) ))                                                     
            ax(iPlot) = gca;
            if numel(stimShading) >= nPlots
                shadeArg = stimShading{iPlot};
            else
                shadeArg = stimShading{1};
            end
            rgbStimShadeColors = [];
            for iEpoch = 1:numel(stimShadingColors)
                rgbStimShadeColors(iEpoch, :) = rgb(stimShadingColors{iEpoch});
            end
            plot_ROI_data(ax(iPlot), currDffAvg, 'EventShading', shadeArg, ...
                                                 'EventShadeColor', rgbStimShadeColors, ...
                                                 'SingleTrials', singleTrials, ...
                                                 'SingleTrialAlpha', singleTrialAlpha, ...
                                                 'StdDevShading', stdDevShading, ...
                                                 'OutlierSD', 4, ...
                                                 'VolumeRate', volumeRate, ...
                                                 'SmoothWinSize', 3);
            yL{iPlot} = ylim(ax(iPlot));
        end
        ax(1).Title.String = plotTitle;

        % Scale y-axes to match
        yMin = min([yL{:}]);
        yMax = max([yL{:}]);
        for iPlot = 1:nPlots
            ylim(ax(iPlot), [yMin yMax]);
        end

        % Save figure -------------------------------------------------------------------------------
        if saveDir
            fileName = [fileNamePrefix, 'ROI_', num2str(iROI), '_', expDate];
            export_fig(fullfile(saveDir, fileName), '-png', f);
            if ~isdir(fullfile(saveDir, 'figFiles'))
                mkdir(fullfile(saveDir, 'figFiles'))
            end
            savefig(f, fullfile(saveDir, 'figFiles', fileName));

        end
    end
catch foldME; rethrow(foldME); end    
    %% PLOT MEAN ROI dF/F THROUGHOUT TRIAL FOR EARLY VS LATE TRIALS
try
    saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
    
    plotTitle = regexprep([expDate, ' - All stims - early vs late trials'], '_', '\\_');
    fileNamePrefix = 'EarlyVsLateTrials_';
    
    s = goodTrials;

    % GROUP BY EARLY/LATE
    groupNames = [];
    groupBounds = [1, 30, 60];
    trialGroups = zeros(1, nTrials);
    for iBound = 1:numel(groupBounds)-1
        trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
        groupNames{iBound} = ['Trials ', num2str(groupBounds(iBound)), '-', num2str(groupBounds(iBound + 1))];
    end
    trialGroups(groupBounds(end):end) = numel(groupBounds);
    groupNames{end + 1} = ['Trials ', num2str(groupBounds(end)), '-', num2str(nTrials)];

    trialGroups = trialGroups(logical(s .* goodTrials));

    singleTrials = 0;
    singleTrialAlpha = 0.3;
    stdDevShading = 1;
    outlierSD = 5;
    legendStr = groupNames;

    ROIlist = 1:size(ROIDffAvg, 3)-1;
%     ROIlist = [1 2 3];

    for iROI = ROIlist       

        % Create figure
        f = figure(iROI); clf
        f.Position = [-1450 50 1000 900];
        f.Color = [1 1 1];

        % Plot reference image and ROI outline
        currPlane = ROImetadata{iROI}(1).plane;
        xData = ROImetadata{iROI}(1).xi;
        yData = ROImetadata{iROI}(1).yi;
        subaxis(2, 3, [1 2], 'ML', 0.06, 'MR', 0.015, 'MT', 0.05, 'MB', 0.08, 'SH', 0)
        hold on
        imshow(infoStruct.refImg{currPlane}, [0 infoStruct.MAX_INTENSITY]);
        plot(xData, yData, 'Color', 'g');
%         title(['Plane #', num2str(infoStruct.ROImetadata{iROI}(1).plane)])

        % Create plot
        currDffAvg = ROIDffAvg(:,logical(s .* goodTrials), iROI); % --> [volume, trial]
        subaxis(2, 3, [4:6])
        ax = gca;
        [~, idx] = max(cellfun(@size, stimShading, repmat({1}, 1, numel(stimShading))));
        plot_ROI_data(ax, currDffAvg, 'EventShading', stimShading{idx}, ...
                                      'TrialGroups', trialGroups,   ...
                                      'SingleTrials', singleTrials, ...
                                      'SingleTrialAlpha', singleTrialAlpha, ...
                                      'OutlierSD', outlierSD, ...
                                      'Legend', legendStr, ...
                                      'VolumeRate', volumeRate, ...
                                      'StdDevShading', stdDevShading);
        title(plotTitle);

        % Save figure -------------------------------------------------------------------------------
        if saveDir
            fileName = [fileNamePrefix, 'ROI_', num2str(iROI), '_', expDate];
            export_fig(fullfile(saveDir, fileName), '-png', f);
            if ~isdir(fullfile(saveDir, 'figFiles'))
                mkdir(fullfile(saveDir, 'figFiles'))
            end
            savefig(f, fullfile(saveDir, 'figFiles', fileName));
        end
    end
catch foldME; rethrow(foldME); end
    %% PLOT MEAN ROI dF/F THROUGHOUT TRIAL COLOR CODED BY BEHAVIOR
try
    saveDir = uigetdir(['D:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
    plotTitle = regexprep([expDate, make_plotTitleSuffix(stimNames)], '(?<!\\)_', '\\_');    
    fileNamePrefix = 'Behavior_Coded_Whole_Trial_Responses_';
    s = infoStruct.stimSepTrials;
    annotValues = [3 0];
    annotTypeInd = 4;
    
    filterVecs = [];
    for iStim = 1:numel(stimGroupNames)
        filterVecs = [filterVecs; s.(stimGroupNames{iStim})];
    end
    trialFilterVecs = logical(filterVecs .* repmat(goodTrials, size(filterVecs, 1), 1));
    
    singleTrials = 1;
    singleTrialAlpha = 0.4;
    stdDevShading = 0;
    
    ROIlist = 1:size(ROIDffAvg, 3)-1;
%     ROIlist = [1 2 3];
    
    for iROI = ROIlist       

        nPlots = size(trialFilterVecs, 1);
        nRows = nPlots + 1;
        
        % Create figure
        f = figure(iROI); clf
        f.Position = [-1450 50 1000 900];
        f.Color = [1 1 1];
        
        % Plot reference image and ROI outline
        currPlane = ROImetadata{iROI}(1).plane;
        xData = ROImetadata{iROI}(1).xi;
        yData = ROImetadata{iROI}(1).yi;
        subaxis(nRows, 3, [1 2], 'ML', 0.06, 'MR', 0.015, 'MT', 0.05, 'MB', 0.08, 'SH', 0)
        hold on
        imshow(infoStruct.refImg{currPlane}, [0 infoStruct.MAX_INTENSITY]);
        plot(xData, yData, 'Color', 'g');
%         title(['Plane #', num2str(infoStruct.ROImetadata{iROI}(1).plane)])
        
        clear ax
        yL = [];
        for iPlot = 1:nPlots
            currDffAvg = ROIDffAvg(:, trialFilterVecs(iPlot, :), iROI); % --> [volume, trial]
            annotData = annotationTypes{annotTypeInd}.volAnnotArr;
            currAnnotArr = annotData(trialFilterVecs(iPlot, :), :);
            
            subaxis(nRows, 3, ( (iPlot*3 + 1):(iPlot*3 + 3) ))
            ax(iPlot) = gca;
            if numel(stimShading) >= nPlots
                shadeArg = stimShading{iPlot};
            else
                shadeArg = stimShading{1};
            end
            rgbStimShadeColors = [];
            for iEpoch = 1:numel(stimShadingColors)
                rgbStimShadeColors(iEpoch, :) = rgb(stimShadingColors{iEpoch});
            end
            plot_ROI_data(ax(iPlot), currDffAvg, 'AnnotArray', currAnnotArr', ...
                'AnnotValues', annotValues', ...
                'EventShading', shadeArg, ...
                'EventShadeColor', rgbStimShadeColors, ...
                'SingleTrials', singleTrials, ...
                'SingleTrialAlpha', singleTrialAlpha, ...
                'StdDevShading', stdDevShading, ...
                                                 'VolumeRate', volumeRate, ...
                                                 'OutlierSD', 4);            
            yL{iPlot} = ylim(ax(iPlot));
        end
        ax(1).Title.String = plotTitle;

        % Scale y-axes to match
        yMin = min([yL{:}]);
        yMax = max([yL{:}]);
        for iPlot = 1:nPlots
            ylim(ax(iPlot), [yMin yMax]);
        end

        % Save figure -------------------------------------------------------------------------------
        if saveDir
            fileName = [fileNamePrefix, 'ROI_', num2str(iROI), '_Plane_', num2str(currPlane), '_', expDate];
            export_fig(fullfile(saveDir, fileName), '-png', f);
            if ~isdir(fullfile(saveDir, 'figFiles'))
                mkdir(fullfile(saveDir, 'figFiles'))
            end
            savefig(f, fullfile(saveDir, 'figFiles', fileName));
            
        end
    end
catch foldME; rethrow(foldME); end    
%%
try
% Average and downsample ficTrac moveSpeed
smMoveSpeed = movmean(movmean(moveSpeed * FRAME_RATE, 11, 1), 11, 1);
volMoveSpeed = smMoveSpeed(volFrames, :);
ROIdata = ROIDffAvg(:,:,1);

cm = parula(numel(unique(volMoveSpeed(:))));
msSorted = sort(unique(volMoveSpeed(:)));


figure(1);clf;hold on

for iTrial = 1:nTrials
    for iVol = 1:nVolumes-1
        currVal = volMoveSpeed(iVol, iTrial);
        [~, test] = min(abs(msSorted - currVal));
        plot([iVol, iVol + 1], [ROIdata(iVol, iTrial), ROIdata(iVol + 1, iTrial)], 'color', cm(test,:)); 
    end
end
catch foldME; rethrow(foldME); end

%% %=================================================================================================
%%% FICTRAC DATA PLOTTING
%% %================================================================================================
s = stimSepTrials;

currStim = logical(stimTrialGroups);
smWin = 3;
nBins = 50;
thresh = 0.3;
nROIs = size(ROIDataAvg, 3)-3;
currTrial = 97;

try
    
% Pull out good trials from data 
mmSpeedData = ftData.moveSpeed(:,logical(goodTrials .* currStim)) * FRAME_RATE * 4.5;   % --> [frame, trial] (mm/sec)
dHD = rad2deg(ftData.yawSpeed(:,logical(goodTrials .* currStim)) * FRAME_RATE);         % --> [frame, trial] (deg/sec)
dHDSmooth = movmean(movmean((dHD), smWin, 1), smWin, 1);                                 
goodFl = ROIDffAvg(:,logical(goodTrials .* currStim),:);
% goodFl = ROIDataAvg(:,logical(goodTrials .* currStim),:);   % --> [volume, trial, ROI]

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
allSpeed = as_vector(movmean(movmean(mmSpeedData(volFrames,:), smWin, 2), smWin, 2));
allYaw =  as_vector(dHDSmooth(volFrames,:)); % Note that this is currently directionless yaw speed

corrMat = [allSpeed, allYaw];
for iROI = 1:nROIs
    corrMat(:, end + 1) = as_vector(movmean(goodFl(:,:,iROI), smWin, 2)); % --> columns: [speed, yaw, ROI1, ROI2, ...]
end

% % Find best fitting # of lag frames for each ROI
% [xCorrROI1, ~] = xcorr(as_vector(mmSpeedData(volFrames,:)), as_vector(goodFl(:,:,1)), 'coeff');
% [~, indROI1] = max(xCorrROI1);
% [xCorrROI2, lags] = xcorr(as_vector(mmSpeedData(volFrames,:)), as_vector(goodFl(:,:,2)), 'coeff');
% [~, indROI2] = max(xCorrROI2);
% bestLags = lags([indROI1, indROI2]);
% % bestLags(bestLags < 0) = 0;
bestLags = [0 0];

% Calculate corrCoeffs for each ROI
allTrialCorrCoeffs = [];
for iROI = 1:nROIs
    currCorrMat = corrcoef(corrMat(1:end, 1), corrMat(1:end, iROI + 2));
    allTrialCorrCoeffs(iROI) = currCorrMat(1,2);
end
disp(allTrialCorrCoeffs);

% ---------------------------------------------------------------------------------------------------
%  Full dataset plots
% ---------------------------------------------------------------------------------------------------

plotNum = 1;
if ~exist('figHandles')
    figHandles = [];
end

% Remove frames if move speed falls below threshold
dffDataVec = []; flDataThresh = [];
mmSpeedDataSmooth = movmean(movmean(mmSpeedData(volFrames, :), smWin, 1), smWin, 1);
flDataSmooth = movmean(movmean(goodFl, smWin, 1), smWin, 1);
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
    title(['ROI ', num2str(iROI), '  (r = ', sprintf('%.3f)', allTrialCorrCoeffs(iROI))]); xlabel('Move speed (mm/sec)'); ylabel('Raw F (AU)')
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
histogram2(allYaw, smooth(as_vector(mmSpeedData(volFrames, :))), nBins, 'DisplayStyle', 'tile')
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
imagesc(movmean(speedDataTrim', smWin, 2));
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

shadeVols = [];
if ~isempty(infoStruct.stimOnsetTimes)
    shadeTimes = [infoStruct.stimOnsetTimes(1), infoStruct.stimOnsetTimes(1) + infoStruct.stimDurs(1)];
    shadeVols = round(shadeTimes * volumeRate);
end
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
currMoveSpeed = mmSpeedData(volFrames,currTrial);
yyaxis left
plot(smooth(smooth(currMoveSpeed,smWin),smWin), 'LineWidth', 2, 'Color', plotColors(1,:))
ylabel('Move speed (mm/sec)');
lgdStr = {'Move speed'};
ylim([0 max(mmSpeedData(:))]);
% ylim([0 30])
yyaxis right
for iROI = 1:nROIs
    plot(smooth(goodFl(:,currTrial,iROI),smWin), ':', 'LineWidth', 2, 'Color', plotColors(iROI + 1, :))
    lgdStr{end + 1} = ['ROI ', num2str(iROI)];
end
ylabel('Raw florescence (AU)');
legend(lgdStr, 'autoupdate', 'off');
if ~isempty(infoStruct.stimOnsetTimes)
    plot_stim_shading(shadeVols, 'Axes', ax);
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
ylim([0 max(goodFl(:))])
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
plot(smooth(smooth(currHD,smWin),smWin), 'LineWidth', 2, 'Color', plotColors(1,:))
ylabel('Yaw speed (deg/sec)');
lgdStr = {'Yaw speed'};
yyaxis right
for iROI = 1:nROIs
    plot(smooth(goodFl(:,currTrial,iROI),smWin), '--', 'LineWidth', 2, 'Color', plotColors(iROI + 1, :))
    lgdStr{end + 1} = ['ROI ', num2str(iROI)];
end
ylabel('Raw florescence (AU)');
legend(lgdStr, 'autoupdate', 'off');
if ~isempty(infoStruct.stimOnsetTimes)
    plot_stim_shading(shadeVols, 'Axes', ax);
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
    plot(smooth(currMoveSpeed,smWin), smooth(goodFl(:, currTrial, iROI),smWin), 'o', 'color', 'b', 'markersize', 5);
    C = corrcoef(smooth(currMoveSpeed, smWin), smooth(goodFl(:, currTrial, iROI), smWin));
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
    smoothVelData = smooth(smooth(currMoveSpeed, smoothWin), smoothWin);  % --> [sample, axis, trial]

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
        plot(axMove, smooth(currIntXY(currSegInd, 1), 3), smooth(currIntXY(currSegInd, 2), 3), 'Color', cm(iSeg, :), 'LineWidth', 2)
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
    plot(smooth(currYawSpeed, 11), 'LineWidth', 2, 'color', 'r');
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
planeNum = 12;
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
        plot(smooth(concatROIDff(:,iROI), 3), 'b');
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
    currX = smooth(mmXYdata(:,1,iTrial), smWin);
    currY = smooth(mmXYdata(:,2,iTrial), smWin);
    currSpeed = smooth(mmSpeedData(:,iTrial), smWin);
    currFWSpeed = smooth(mmFWSpeed(:,iTrial), smWin);
    currYawSpeed = smooth(dHD(:,iTrial), smWin);
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
        fwSpeedSmooth = smooth(smooth(currFWSpeed(:), smWin), smWin);
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
            currFl = smooth(goodFlNorm(:,iTrial,iROI),smWin);
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
        yawSpeedSmooth = smooth(smooth(currYawSpeed,smWin),smWin);
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
            currFl = smooth(goodFlNorm(:,iTrial,iROI),smWin);
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