%% LOAD REGISTERED DATA FILE AND BEHAVIORAL ANNOTATION DATA

[dataFile, parentDir, ~] = uigetfile('*.mat', 'Select an imaging data session file', 'B:\Dropbox (HMS)\2P Data\Imaging Data\');

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
        
        % Load reference images file
        disp('Loading reference images...')
%         [refImgFile, refImgPath, ~] = uigetfile('*.mat', 'Select a reference image file', parentDir);
%         refImgPath = parentDir; refImgFile = ['sid_', num2str(analysisMetadata.sid), '_refImages.mat'];
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
        disp('Loading FicTrac data...')
        ftDir = fullfile('B:\Dropbox (HMS)\2P Data\Behavior Vids\', analysisMetadata.expDate, '_Movies\FicTracData');
        if isdir(ftDir)
            ftData = load_fictrac_data('Sid', analysisMetadata.sid, 'AnalysisMetadata', fullfile(parentDir, 'analysisMetadata.mat'), 'ParentDir', ftDir);
        else
            ftData = load_fictrac_data('Sid', analysisMetadata.sid, 'AnalysisMetadata', fullfile(parentDir, 'analysisMetadata.mat'));
        end
        if ~isfield(analysisMetadata, 'ftData')
            analysisMetadata.ftData = ftData;
            save(fullfile(parentDir, 'analysisMetadata.mat'), 'analysisMetadata', '-v7.3');
        end
        
        % Load volume-averaged raw session data
        volAvgDataFile = fullfile(parentDir, ['sid_', num2str(analysisMetadata.sid), '_volAvgSessionData.mat']);
        if exist(volAvgDataFile)
            load(volAvgDataFile) % "volAvgSessionData"
        end
        
        % Reload analysis metadata in case .goodTrials has been updated
        load(fullfile(parentDir, 'analysisMetadata.mat')) % --> 'analysisMetadata' info struct
        analysisMetadata.refImg = refImages;
        disp('All data loaded')
        
        % Omit any trials in which FicTrac reset from analysis
        analysisMetadata.goodTrials(logical(ftData.resets)) = 0;
        
        % ------- Copy variables for convenience -------
        sessionSize = size(m, 'wholeSession');
        expDate = analysisMetadata.expDate;
        sid = analysisMetadata.sid;
        nPlanes = analysisMetadata.nPlanes;
        nVolumes = analysisMetadata.nVolumes;
        refImg = analysisMetadata.refImg;
        if ~isempty(analysisMetadata.nFrames)
            nFrames = analysisMetadata.nFrames;
        else
            nFrames = nVolumes;
        end
        nTrials = analysisMetadata.nTrials;
        nGoodTrials = sum(analysisMetadata.goodTrials);
        stimTypes = analysisMetadata.stimTypes;
        stimOnsetTimes = analysisMetadata.stimOnsetTimes;
        stimDurs = analysisMetadata.stimDurs;
        trialDuration = analysisMetadata.trialDuration;
        volumeRate = analysisMetadata.volumeRate;
        volFrames = analysisMetadata.volFrames;
        goodTrials = analysisMetadata.goodTrials;
        stimSepTrials = analysisMetadata.stimSepTrials;
        behaviorAnnotArr = annotationTypes{contains(annotationTypeSummary.AnnotationType, 'move')}.frameAnnotArr;
        
        % Create hardcoded parameters
        FRAME_RATE = 25; % This is the frame rate for the behavior video
        MAX_INTENSITY = analysisMetadata.MAX_INTENSITY;
        if isempty(nFrames)
            nFrames = sum(trialDuration) * FRAME_RATE;
        end
        volTimes = (1:nVolumes)' ./ volumeRate;
        frameTimes = (1:nFrames)' ./ FRAME_RATE;
        
        % Create directory for saving analysis files if necessary
        if ~isdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'])
            mkdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'])
        end
        
    end%if
catch(foldME); rethrow(foldME); end
%% Plot raw fluorescence minimum for each trial
try
rawFlAvg = squeeze(min(min(volAvgSessionData, [], 1), [], 2)); % --> [plane, trial]
% rawFlAvg = squeeze(mean(mean(volAvgSessionData, 1), 2)); % --> [plane, trial]


flThresh = []; 
omitTrials = [];
% omitTrials = [21 44 45 61 98 99 110];
% omitTrials = [8 9 10 12 13 14 17 18 19 23 30:39 55 58 138];
find(rawFlAvg(5,:) > 150 | rawFlAvg(5,:) < 50);

if ~isempty(omitTrials)
    rawFlAvg(:,omitTrials) = [];
end
figure(1);clf;
plot(rawFlAvg');
if ~isempty(flThresh)
   ylim([0 flThresh]) 
end

catch foldME; rethrow(foldME); end
%% Remove trials from analysis
try
    
goodTrials = analysisMetadata.goodTrials;
goodTrials(omitTrials) = 0;

catch(foldME); rethrow(foldME); end

%% =================================================================================================
%   BEHAVIOR SUMMARIES                                   
%%%=================================================================================================
    %% PLOT 2-D SUMMARY OF BEHAVIOR DATA ANNOTATIONS
try
saveFig = uigetdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
s = stimSepTrials;

stimNames = {'EtOH\_neat', 'EtOH\_e-2'};

% 
trialGroups = [];
plotTitleSuffix = '';
fileNameSuffix = '_AllTrials';
plotAnnotTypes = [2]; % 1 = stims, 2 = behavior
% 
% trialGroups = [s.OdorA + 2 * s.OdorB]; 
% plotTitleSuffix = ['  —  ', stimNames{1}, ' (top) vs. ', stimNames{2}, ' (bottom)']; %
% fileNameSuffix = ['_OdorAvsOdorB'];
% plotAnnotTypes = [2]; % 1 = odor stims, 2 = behavior

try

% Create plot titles
nPlots = length(plotAnnotTypes);
titleStrings = [];
plotNames = [];
annotArr = [];

for iPlot = 1:nPlots
    if plotAnnotTypes(iPlot) == 1
        
        stimAnnotTypes = annotationTypes(contains(annotationTypeSummary.AnnotationType, analysisMetadata.stimTypes));        
        tempAnnotArr = zeros(size(stimAnnotTypes{1}.frameAnnotArr));
        for iType = 1:numel(stimAnnotTypes)
           tempAnnotArr = tempAnnotArr + (iType * stimAnnotTypes{iType}.frameAnnotArr); 
        end
        tempAnnotArr = tempAnnotArr > 0;
        
        plotNames{iPlot} = 'Stim Delivery';
        titleStrings{iPlot} = [regexprep(expDate, '_', '\\_'), '    ', [plotNames{iPlot}, ' summary ', plotTitleSuffix]]; % regex to add escape characters
        annotArr{iPlot} = tempAnnotArr;
    elseif plotAnnotTypes(iPlot) == 2
        % Behavior
        plotNames{iPlot} = 'Behavior Annotation';
        titleStrings{iPlot} = [regexprep(expDate, '_', '\\_'), '    ', [plotNames{iPlot}, ' summary ', plotTitleSuffix]]; % regex to add escape characters
        annotArr{iPlot} = behaviorAnnotArr;
    end%if
end%for

% Create figure
f = figure(7);clf
if nPlots == 1
    f.Position = [-1050 45 1020 950];
elseif nPlots == 2
    f.Position = [-1050 45 900 950];
elseif nPlots == 3
    f.Position = [-1050 45 700 950];
end
f.Color = [1 1 1];

% Create and format each plot
for iPlot = 1:nPlots
    
    subaxis(nPlots, 1, iPlot, ...
            'MarginTop', 0, ...
            'MarginBottom', 0.055, ...
            'MarginRight', 0.015, ...
            'MarginLeft', 0.065, ...
            'Spacing', 0, ...
            'PaddingTop', 0.03 ...
            );
    ax = gca;
    [~, ax, ~] = plot_behavior_summary_2D(analysisMetadata, annotArr{iPlot}, ax, titleStrings{iPlot}, trialGroups);
    ax.FontSize = 14;
    ax.Title.FontSize = 12;
    ax.XLabel.FontSize = 14;
    if iPlot ~= nPlots
        ax.XLabel = [];
        ax.XTickLabels = [];
    end 
    
    % Plot stim times
    hold on
    stimOnsetFrame = stimOnsetTimes(1) * FRAME_RATE;
    stimOffsetFrame = (stimOnsetTimes(1) + stimDurs(1)) * FRAME_RATE;
    plot(ax, [stimOnsetFrame, stimOnsetFrame], ylim(), 'Color', 'k', 'LineWidth', 2)
    plot(ax, [stimOffsetFrame, stimOffsetFrame], ylim(), 'Color', 'k', 'LineWidth', 2)
end

if saveFig
    % Create analysis directory if necessary
    saveDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'];
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    
    % Create filename
%     for iPlot = 1:nPlots
%         fileName = [regexprep(plotNames{iPlot}, ' ', '')];
%     end
    fileName = ['Behavior_Annotation_Summary ', fileNameSuffix, '_', expDate];
    
    % Warn user and offer to cancel save if this will overwrite existing files
    overwrite = 1;
    if exist(fullfile(saveDir, [fileName, '.fig']), 'file') ~= 0 || ...
            exist(fullfile(saveDir, [fileName, '.png']), 'file') ~= 0 
        dlgAns = questdlg('Saving this figure will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
        if strcmp(dlgAns, 'No')
            overwrite = 0;
            disp('Saving cancelled')
        end
    end
    
    % Save figure files
    if overwrite
        export_fig(fullfile(saveDir, fileName), '-png', f);
        if ~isdir(fullfile(saveDir, 'figFiles'))
            mkdir(fullfile(saveDir, 'figFiles'))
        end
        savefig(f, fullfile(saveDir, 'figFiles', fileName));
    end
end%if

clear saveFig plotTypes trialGroups plotTitleSuffix fileNameSuffix s nPlots titleStrings plotNames annotArr f saveDir fileName overwrite dlgAns
catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end
    %% PLOT 1-D SUMMARY OF BEHAVIOR DATA ANNOTATIONS
try
%----- Plot 1D trial-averaged movement data -----
s = stimSepTrials; 
actionNames = {'NA', 'Locomotion', 'Grooming', 'IsoMovement'};
% stimNames = {'Ethanol\_e-1', 'Ethanol\_e-3'};

saveFig = uigetdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
actionNum = [2]; % locomotionLabel = 2; noActionLabel = 0; groomingLabel = 3; isoMovementLabel = 4;
actionName = actionNames{actionNum};
figTitle = regexprep([expDate, '  —  Fly ', actionName, ' throughout trial (red = stim period)'], '_', '\\_');
% 
% ALL TRIALS
trialGroups = [goodTrials];
fileNameSuffix = ['_AllTrials_', actionName];
groupNames = {'All trials'};

% % GROUP BY STIM TYPE
% trialGroups = [s.OdorA + 2 * s.OdorB] .* goodTrials; 
% fileNameSuffix = ['_OdorAvsOdorB_', actionName]; 
% groupNames = stimNames;

% % GROUP BY EARLY/LATE
% groupNames = [];
% groupBounds = [1, 35];
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
%    groupNames{iBound} = ['Trials ', num2str(groupBounds(iBound)), '-', num2str(groupBounds(iBound + 1))];
% end
% trialGroups(groupBounds(end):end) = numel(groupBounds);
% groupNames{end + 1} = ['Trials ', num2str(groupBounds(end)), '-', num2str(nTrials)];
% fileNameSuffix = '_EarlyVsLateTrials';
% % 

stimShadingColors = {'red', 'green'};

try
   
% Get odor stim times
odorOnset = mode(analysisMetadata.stimOnsetTimes(logical(s.OdorA + s.OdorB)));
odorOffset = odorOnset + mode(analysisMetadata.stimDurs(logical(s.OdorA + s.OdorB)));
odorTimes = [odorOnset, odorOffset];
odorFrames = floor(odorTimes * FRAME_RATE);

stimShading = {odorFrames};

% Create array of annotation data
f = figure(2); clf; hold on
f.Position = [100 100 1600 500];
f.Color = [1 1 1];

if isempty(trialGroups)
    
    % Plot summed movement data
    annotArrSum = sum(ismember(behaviorAnnotArr, actionNum), 1) ./ nTrials;
    ax = gca(); 
    ax.FontSize = 14;
    plot_behavior_summary_1D(analysisMetadata, annotArrSum(2:end-1), ax, figTitle);
    
    % Add shading during stimulus presentations
    yL = ylim();
    for iType = 1:numel(stimShading)
        for iStim = 1:size(stimShading{iType}, 1)
            stimStart = stimShading{iType}(iStim, 1);
            stimLength = stimShading{iType}(iStim, 2) - stimShading{iType}(iStim, 1);
            rectPos = [stimStart, yL(1), stimLength, diff(yL)]; % [x y width height]
            rectangle('Position', rectPos, 'FaceColor', [rgb(stimShadingColors{iType}), 0.1], 'EdgeColor', 'none');
            ylim(yL);
        end
    end
else
    annotArrSum = [];
    yLimsAll = [];
    ax = [];
    nGroups = length(unique(trialGroups(trialGroups ~= 0)));
    cm = parula(nGroups);
    cm = [rgb('blue'); rgb('green'); rgb('red'); rgb('magenta'); rgb('cyan')];
    for iGroup = 1:nGroups
        
        % Plot summed movement data
%         if length(unique(trialGroups(trialGroups ~= 0))) > 1
%             f.Position = [100 50 1000 950];
%         end
%         ax{iGroup} = subplot(nGroups, 1, iGroup);
        ax = gca();
        ax.FontSize = 14;
        colormap(jet(nGroups))
        annotArrSum = sum(ismember(behaviorAnnotArr(trialGroups == iGroup, :), actionNum), 1) ./ sum(trialGroups == iGroup);
        [plt, ~, ~] = plot_behavior_summary_1D(analysisMetadata, annotArrSum, 'PlotAxes', ax, 'LineColor', cm(iGroup, :));
        plt.LineWidth = 2;
        if iGroup ~= length(unique(trialGroups))
            xlabel('');
        end
        
        % Add shading during stimulus presentations
        yL = ylim();
        yLimsAll(iGroup, :) = yL;
        for iType = 1:numel(stimShading)
            for iStim = 1:size(stimShading{iType}, 1)
                stimStart = stimShading{iType}(iStim, 1);
                stimLength = stimShading{iType}(iStim, 2) - stimShading{iType}(iStim, 1);
                rectPos = [stimStart, 0, stimLength, 1000]; % using large height value in case yLims increase later
                rectangle('Position', rectPos, 'FaceColor', [rgb(stimShadingColors{iType}), 0.05], 'EdgeColor', 'none');
                ylim(yL);
            end
        end
    end%iGroup
    legend(groupNames, 'FontSize', 14, 'Location', 'Best')
    ax.XLim = [20 nFrames-20]; % to improve plot appearance
    ax.YLim = [0 1];
%     % Make sure all plots use the same yLims
%     yLimMax = max(yLimsAll(:));
%     for iGroup = 1:length(unique(trialGroups(trialGroups~=0)))
%         ylim(ax{iGroup}, [yLimsAll(iGroup, 1), yLimMax]);
%     end
    suptitle(figTitle);
end

if saveFig
    % Create analysis directory if necessary
    saveDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'];
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    
    % Create filename
    fileName = regexprep(['Summed_Movement', fileNameSuffix, '_', expDate], '_', '\_');
    
    % Warn user and offer to cancel save if this will overwrite existing files
    overwrite = 1;
    if exist(fullfile(saveDir, [fileName, '.fig']), 'file') ~= 0 || ...
            exist(fullfile(saveDir, [fileName, '.png']), 'file') ~= 0 
        dlgAns = questdlg('Saving this figure will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
        if strcmp(dlgAns, 'No')
            overwrite = 0;
            disp('Saving cancelled')
        end
    end
    
    % Save figure files
    if overwrite
        export_fig(fullfile(saveDir, fileName), '-png', f);
        if ~isdir(fullfile(saveDir, 'figFiles'))
            mkdir(fullfile(saveDir, 'figFiles'))
        end
        savefig(f, fullfile(saveDir, 'figFiles', fileName));
    end
end%if

clear s saveFig fileNameSuffix actionLabel trialGroups figTitle plotNames stimShadingColors odorOnset odorOffset odorTimes odorFrames laserOnset laserOffsetlaserTimes
clear laserFrames stimShading f annotArrSum ax yL stimStart stimLength rectPos yLimsAll yLimMax saveDir fileName overwrite dlgAns
catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end
    %% PLOT 2-D SUMMARY OF FICTRAC DATA
try
saveFig = uigetdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
s = stimSepTrials;

ftVarName = 'moveSpeed'; % 'moveSpeed', 'fwSpeed', 'yawSpeed'
sdCap = 1.5;
smWin = 9;
cmName = @parula;
figTitle = [regexprep(expDate, '_', '\\_'), '  —  FicTrac ', ftVarName];

% % 
% ALL TRIALS
trialGroups = [];
fileNameSuffix = ['_AllTrials'];
figTitleSuffix = '';
% 
% 
% % GROUP BY STIM TYPE
% trialGroups = [[s.OdorA + 2 * s.OdorB] .* goodTrials]; 
% fileNameSuffix = ['_OdorAvsOdorB']; 
% figTitleSuffix = ['  —  ', stimNames{1}, ' (top) vs. ', stimNames{2}, ' (bottom)']

try

% Extract FicTrac data
rawData = analysisMetadata.ftData.(ftVarName);          % --> [frame, trial]
rawData = rawData';                                     % --> [trial, frame]
if strcmp(ftVarName, 'yawSpeed')
    plotData = abs(rad2deg(rawData .* FRAME_RATE));    	% --> [trial, frame] (deg/sec)
    figTitle = [figTitle, ' (deg/sec)'];
else
    plotData = rawData .* FRAME_RATE .* 4.5;            % --> [trial, frame] (mm/sec)
    figTitle = [figTitle, ' (mm/sec)'];
end

% Cap values at 3 SD above mean
capVal = mean(plotData(:)) + (sdCap * std(plotData(:)));
plotData(plotData > capVal) = capVal;

% Smooth data
smPlotData = movmean(plotData, smWin, 2);

% Create colormap
cm = cmName(numel(unique(smPlotData)));
if ~isempty(trialGroups)
    cm = [0 0 0; cm(2:end, :)];
end

% Plot data
titleStr = [figTitle, figTitleSuffix];
fileName = ['FicTrac_', ftVarName '_Summary', fileNameSuffix, '_', expDate];
[~, ax, ~] = plot_2D_summary(analysisMetadata, smPlotData, ...
                'trialGroups', trialGroups, ...
                'titleStr', titleStr, ...
                'saveDir', saveFig, ...
                'fileName', fileName, ...
                'sampRate', FRAME_RATE, ...
                'colormap', cm ...
                );

% Plot stim times
colorbar
hold on
stimOnsetFrame = stimOnsetTimes(1) * FRAME_RATE;
stimOffsetFrame = (stimOnsetTimes(1) + stimDurs(1)) * FRAME_RATE;
plot(ax, [stimOnsetFrame, stimOnsetFrame], ylim(), 'Color', 'k', 'LineWidth', 2)
plot(ax, [stimOffsetFrame, stimOffsetFrame], ylim(), 'Color', 'k', 'LineWidth', 2)
    
catch foldME; rethrow(foldME); end
catch foldME; rethrow(foldME); end
    %% PLOT 1-D SUMMARY OF FICTRAC DATA
try
    
s = stimSepTrials; 
% stimNames = {'Ethanol\_e-1', 'Ethanol\_e-3', '200Hz tone', 'No Stim'};
saveFig = uigetdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');

includeQuiescence = 0;
if ~includeQuiescence
    fileNameSuffix = 'NoQuiescence_';
else
    fileNameSuffix = '';
end
figTitle = [expDate, '  —  Trial-Averaged FicTrac data'];
trialGroups = goodTrials';
smWin = 11;

% % 
% ALL TRIALS
trialGroups = [goodTrials];
fileNameSuffix = [fileNameSuffix, 'AllTrials'];
groupNames = {'All trials'};

% % 
% % GROUP BY STIM TYPE
% trialGroups = [[s.OdorA + 2 * s.OdorB] .* goodTrials]; 
% fileNameSuffix = [fileNameSuffix, 'OdorAvsOdorB']; 
% groupNames = stimNames;
% % % % 
% 
% % GROUP BY EARLY/LATE
% groupNames = [];
% groupBounds = [1, 35];
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
%    groupNames{iBound} = ['Trials ', num2str(groupBounds(iBound)), '-', num2str(groupBounds(iBound + 1))];
% end
% trialGroups(groupBounds(end):end) = numel(groupBounds);
% groupNames{end + 1} = ['Trials ', num2str(groupBounds(end)), '-', num2str(nTrials)];
% fileNameSuffix = [fileNameSuffix, 'EarlyVsLateTrials'];
% % 

% % GROUP BY EARLY/LATE FOR A SINGLE STIM TYPE
% targetStim = s.OdorA;
% stimName = 'OdorA';
% groupNames = [];
% groupBounds = [1, 30, 70];
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
%    groupNames{iBound} = [stimName, ' trials ', num2str(groupBounds(iBound)), '-', num2str(groupBounds(iBound + 1))];
% end
% trialGroups(groupBounds(end):end) = numel(groupBounds);
% trialGroups(logical(~targetStim .* goodTrials)) = 0;
% groupNames{end + 1} = [stimName, ' trials ', num2str(groupBounds(end)), '-', num2str(nTrials)];
% fileNameSuffix = [fileNameSuffix, 'EarlyVsLate_', stimName];

try
    
% Extract relevant data
shadeFrames = [];
if ~isempty(analysisMetadata.stimOnsetTimes)
    shadeTimes = [analysisMetadata.stimOnsetTimes(1), analysisMetadata.stimOnsetTimes(1) + analysisMetadata.stimDurs(1)];
    shadeFrames = round(shadeTimes * FRAME_RATE);
end
xTickFR = [0:1/trialDuration:1] * (trialDuration * FRAME_RATE);
xTickLabels = [0:1/trialDuration:1] * trialDuration;
mmSpeedData = ftData.moveSpeed * FRAME_RATE * 4.5;  % --> [frame, trial] (mm/sec)
dHD = abs(rad2deg(ftData.yawSpeed * FRAME_RATE));        % --> [frame, trial] (deg/sec)
fwSpeed = ftData.fwSpeed * FRAME_RATE * 4.5;        % --> [frame, trial  (mm/sec)
nFrames = size(mmSpeedData, 1);
stimShadingColors = {'red', 'green'};

% Create figure
f = figure(3); clf; hold on
f.Position = [-1100 50 900 930];
f.Color = [1 1 1];

% Create axes
M = 0.02;
P = 0.00;
axVel = subaxis(3,1,1, 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.05); hold on
axFWSpeed = subaxis(3,1,2, 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.06); hold on
axYawSpeed = subaxis(3,1,3, 'S', 0, 'M', M, 'PB', 0.06, 'PL', 0.06); hold on

% Plot data
nGroups = length(unique(trialGroups(trialGroups ~= 0)));
cm = [rgb('blue'); rgb('red'); rgb('green'); rgb('magenta'); rgb('cyan')];
for iGroup = 1:nGroups
    
    % Calculate mean values for current group
    currXYSpeed = mmSpeedData(:, trialGroups==iGroup);
    currFWSpeed = fwSpeed(:, trialGroups == iGroup);
    currYawSpeed = dHD(:, trialGroups == iGroup);
    
    if ~includeQuiescence
        currAnnotData = behaviorAnnotArr';
        currAnnotData(:, trialGroups ~= iGroup) = [];
        currXYSpeed(currAnnotData ~= 2) = nan;
        currFWSpeed(currAnnotData ~= 2) = nan;
        currYawSpeed(currAnnotData ~= 2) = nan;
    end
    
    % Omit outliers
    outlierCalc = @(x) mean(x) + 4 * std(x);
    currXYSpeed(currXYSpeed >= outlierCalc(mmSpeedData(:))) = nan;
    currFWSpeed(currFWSpeed >= outlierCalc(fwSpeed(:))) = nan;
    currYawSpeed(currYawSpeed >= outlierCalc(dHD(:))) = nan;
%     
    meanSpeed = smooth(mean(currXYSpeed, 2, 'omitnan'), smWin);
    meanFWSpeed = smooth(mean(currFWSpeed, 2, 'omitnan'), smWin);
    meanYawSpeed = smooth(mean(currYawSpeed, 2, 'omitnan'), smWin);
    
%     % Calculate mean values for current group
%     currXYSpeed = mmSpeedData(:, trialGroups==iGroup);
%     meanSpeed = smooth(mean(currXYSpeed, 2), smWin);
%     meanFWSpeed = smooth(mean(fwSpeed(:,trialGroups==iGroup), 2), smWin);
%     meanYawSpeed = smooth(mean(dHD(:,trialGroups==iGroup), 2), smWin);
    
    % XY speed plot
    axes(axVel)
    plot(meanSpeed, 'linewidth', 2, 'color', cm(iGroup, :));
    
    % Forward speed plot
    axes(axFWSpeed)
    plot(meanFWSpeed, 'linewidth', 2, 'color', cm(iGroup,:));
        
    % Yaw speed plot
    axes(axYawSpeed)
    plot(meanYawSpeed, 'linewidth', 2, 'color', cm(iGroup,:));
    axYawSpeed.XLabel.String = 'Time (sec)';

end%iGroup

% Format axes
axVel.XTick = xTickFR;
axVel.XTickLabel = xTickLabels;
axVel.YLabel.String = 'XY Speed (mm/sec)';
axVel.FontSize = 14;
legend(axVel, groupNames, 'FontSize', 12, 'Location', 'best', 'AutoUpdate', 'off')
axVel.XLim = [9 nFrames-5]; % to improve plot appearance
if ~isempty(shadeFrames)
    plot_stim_shading(shadeFrames, 'Axes', axVel);
end

axFWSpeed.XTick = xTickFR;
axFWSpeed.XTickLabel = xTickLabels;
axFWSpeed.YLabel.String = 'FW Vel (mm/sec)';
axFWSpeed.FontSize = 14;
legend(axFWSpeed, groupNames, 'FontSize', 12, 'Location', 'best', 'AutoUpdate', 'off')
axFWSpeed.XLim = [9 nFrames-5]; % to improve plot appearance
if ~isempty(shadeFrames)
    plot_stim_shading(shadeFrames, 'Axes', axFWSpeed);
end

axYawSpeed.XTick = xTickFR;
axYawSpeed.XTickLabel = xTickLabels;
axYawSpeed.YLabel.String = 'Yaw Speed (deg/sec)';
axYawSpeed.FontSize = 14;
legend(axYawSpeed, groupNames, 'FontSize', 12, 'Location', 'best', 'AutoUpdate', 'off')
axYawSpeed.XLim = [9 nFrames-5]; % to improve plot appearance
if ~isempty(shadeFrames)
    plot_stim_shading(shadeFrames, 'Axes', axYawSpeed);
end
suptitle(regexprep([figTitle, '  —  ', fileNameSuffix], '_', '\\_'));
if saveFig
    % Create analysis directory if necessary
    saveDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'];
    if ~isdir(saveDir)
        mkdir(saveDir);
    end
    
    % Create filename
    fileName = regexprep(['FicTrac_Summary_', fileNameSuffix, '_', expDate], '_', '\_');
    
    % Warn user and offer to cancel save if this will overwrite existing files
    overwrite = 1;
    if exist(fullfile(saveDir, [fileName, '.fig']), 'file') ~= 0 || ...
            exist(fullfile(saveDir, [fileName, '.png']), 'file') ~= 0 
        dlgAns = questdlg('Saving this figure will overwrite one or more existing files in this directory...are you sure you want to do this?', 'Warning', 'Yes', 'No', 'No');
        if strcmp(dlgAns, 'No')
            overwrite = 0;
            disp('Saving cancelled')
        end
    end
    
    % Save figure files
    if overwrite

        export_fig(fullfile(saveDir, fileName), '-png', f);
        if ~isdir(fullfile(saveDir, 'figFiles'))
            mkdir(fullfile(saveDir, 'figFiles'))
        end
        savefig(f, fullfile(saveDir, 'figFiles', fileName));
    end
end%if

catch foldME; rethrow(foldME); end

catch foldME; rethrow(foldME); end
    %% PLOT OVERLAID 2D MOVEMENT DATA AND CALC PRE/POST SINUOSITY
try

% stimNames = {'Ethanol\_neat', 'CO2\_4%', '200Hz tone', 'No Stim'};

startTime = 10;
plotLen = 3;
sinWin = 5;
limScalar = 0.8;

s = stimSepTrials;
shadeTimes = [analysisMetadata.stimOnsetTimes(1), analysisMetadata.stimOnsetTimes(1) + analysisMetadata.stimDurs(1)];
shadeFrames = round(shadeTimes * FRAME_RATE);
xTickFR = [0:1/trialDuration:1] * (trialDuration * FRAME_RATE);
xTickLabels = [0:1/trialDuration:1] * trialDuration;
surfPlot = 0;
mmXY = ftData.intXY * 4.5; % --> [frame, var, trial] (mm)
HD = ftData.intHD;    % --> [frame, trial]
nFrames = size(mmXY, 1);
startFrame = startTime * FRAME_RATE;
endFrame = (startTime + plotLen) * FRAME_RATE;
trialGroups = ones(size(mmXY, 3), 1) .* goodTrials';
cm = [rgb('blue'); rgb('red'); rgb('green'); rgb('magenta'); rgb('cyan')];


% % PLOT ALL TRIALS COLORED BY TIME
% surfPlot = 1;
% cm = jet(endFrame - startFrame + 1);

% % 
% GROUP BY STIM TYPE
trialGroups =  [s.OdorA + 2 * s.OdorB] .* goodTrials;
fileNameSuffix = '_OdorAvsOdorBvsNoStim'; 
groupNames = stimNames;
% 
% % GROUP BY PRECEDING STIM TYPE
% trialGroups = [[0, s.OdorA(1:end-1)] + 2 * [0, s.OdorB(1:end-1)] ...
%                 + 3 * [0, s.NoStim(1:end-1)]] .* goodTrials;
% fileNameSuffix = '_After_OdorAvsOdorBvsNoStim'; 
% groupNames = {['After ', stimNames{1}], ['After ', stimNames{2}], ['After ', stimNames{3}]};

% % GROUP BY EARLY/LATE
% groupNames = [];
% groupBounds = [1, 40, 80];
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
%    groupNames{iBound} = ['Trials ', num2str(groupBounds(iBound)), '-', num2str(groupBounds(iBound + 1))];
% end
% trialGroups(groupBounds(end):end) = numel(groupBounds);
% groupNames{end + 1} = ['Trials ', num2str(groupBounds(end)), '-', num2str(nTrials)];
% fileNameSuffix = '_EarlyVsLateTrials';

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
    legendPlots = [0 0 0 0]; legendObj = [];
    for iGroup = 1:nGroups
        for iTrial = 1:nTrials
            if trialGroups(iTrial) == iGroup
                
                plt = plot(plotX(:,iTrial), plotY(:,iTrial), 'color', cm(iGroup, :), 'linewidth', 1);
                
                % Save one plot line from each group to use in legend
                if ~legendPlots(trialGroups(iTrial))
                    legendObj(trialGroups(iTrial)) = plt;
                    legendPlots(trialGroups(iTrial)) = 1;
                end
            end%if
        end%iTrial
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

% CALCULATE PRE- AND POST- STIM SINUOSITY FOR EACH GROUP

% TODO: remove trials in which the fly isn't moving? Might not really happen in most experiments.

sinuosityPre = []; sinuosityPost = []; groupSinuosities = [];
actualDistPre = []; actualDistPost = [];
for iTrial = 1:size(mmXY, 3)
    
    % Calculate straightness before and after stim
    stimOnsetTime = analysisMetadata.stimOnsetTimes(iTrial);
    stimOnsetFrame = stimOnsetTime * FRAME_RATE;
    currXY = movmean(mmXY(:,:,iTrial), 3, 1);
    euc_dist = @(x1, x2, y1, y2) sqrt((x1 - x2)^2 + (y1 - y2)^2);
    
    % Get all individual distances between frames
    frameDist = [];
    for iFrame = 2:size(currXY, 1)
        startX = currXY(iFrame - 1, 1);
        startY = currXY(iFrame - 1, 2);
        endX = currXY(iFrame, 1);
        endY = currXY(iFrame, 2);
        frameDist(iFrame) = euc_dist(startX, endX, startY, endY);
    end
    
    % Calculate frames bounding the analysis window
    preStartFrame = FRAME_RATE * (stimOnsetTime - sinWin);
    postEndFrame = FRAME_RATE * (stimOnsetTime + sinWin);
    
    % Find pre- and post-stim start/end points
    startXpre = currXY(preStartFrame,1);
    startYpre = currXY(preStartFrame,2);
    endXpre = currXY(stimOnsetFrame, 1);
    endYpre = currXY(stimOnsetFrame, 2);
    startXpost = currXY(stimOnsetFrame,1);
    startYpost = currXY(stimOnsetFrame,2);
    endXpost = currXY(postEndFrame, 1);
    endYpost = currXY(postEndFrame, 2);
    
    % Calculate ideal and actual distance traveled pre- and post-stim
    shortestDistPre = euc_dist(startXpre, endXpre, startYpre, endYpre);
    shortestDistPost = euc_dist(startXpost, endXpost, startYpost, endYpost);
    actualDistPre(iTrial) = sum(frameDist(preStartFrame:stimOnsetFrame));
    actualDistPost(iTrial) = sum(frameDist(stimOnsetFrame:postEndFrame));
    
    % Calculate approximate sinuosities
    sinuosityPre(iTrial) = actualDistPre(iTrial) / shortestDistPre;
    sinuosityPost(iTrial) = actualDistPost(iTrial) / shortestDistPost;
    
end%iTrial

% Make boxplot of sinuosities
import iosr.statistics.*
plotData = []; weightData = [];
figure(4);clf

plotTrials = (trialGroups > 0 & sinuosityPre < 3 & sinuosityPost < 3);

plotData(:,:,1) = tab2box(trialGroups(plotTrials), sinuosityPre(plotTrials));
plotData(:,:,2) = tab2box(trialGroups(plotTrials), sinuosityPost(plotTrials));

weightData(:,:,1) = tab2box(trialGroups(plotTrials), actualDistPre(plotTrials));
weightData(:,:,2) = tab2box(trialGroups(plotTrials), actualDistPost(plotTrials));
weightData = ones(size(weightData));
bp = boxPlot(groupNames, plotData, ...
            'weights', weightData, ...
            'showOutliers', true, ...
            'GroupLabels', {'Pre', 'Post'}, ...
            'style', 'hierarchy', ...
            'showscatter', true, ...
            'showviolin', false, ...
            'linecolor', 'k' ...
            , 'showmean', true ...
            , 'notch', true ...
            );

for iGroup = 1:nGroups
    % Divide into trial groups, dropping unrealistically large values, and calculate means weighted
    % by the actual distance traveled during the time window
    preMean = wmean(sinuosityPre(trialGroups==iGroup & sinuosityPre < 100), actualDistPre(trialGroups==iGroup & sinuosityPre < 100));
    postMean = wmean(sinuosityPost(trialGroups == iGroup & sinuosityPost < 100), actualDistPost(trialGroups==iGroup & sinuosityPost < 100));
    groupMeanSinuosity(iGroup, [1 2]) = [preMean, postMean]; 
end

meanSinuosities = table(groupMeanSinuosity(:,1), groupMeanSinuosity(:,2), 'RowNames', groupNames, 'VariableNames', {'BeforeOnset', 'AfterOnset'});
disp(meanSinuosities)

catch(foldME); rethrow(foldME); end
catch(foldME); rethrow(foldME); end

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
catch(foldME); rethrow(foldME); end
%% =================================================================================================
%   ODOR/SOUND STIM HEATMAPS                                
%%%=================================================================================================
try
    % Show summary again
    odorEventName = 'OdorA';
    eventInd = ~cellfun(@isempty, strfind(primaryEventNames, odorEventName));
    currSummary = allCondSummaries{eventInd};
    currCondNames = allCondNames{eventInd};
    disp(currSummary)
    
    % Calculate absolute max dF/F value across all planes and stim types
    currConds = [3 6];
    makeVid = 1;
    sigma = [0.6];   
    rangeType = 'Max';
    rangeScalar = 0.7;
    saveDir = [];
    fileName = [odorEventName, '_Response_Heatmaps'];

    plotTitles = [];
    for iCond = currConds
        plotTitles{iCond} = [regexprep(currCondNames{iCond}, '_', '\\_'), ' ', currSummary{iCond, 3}{:}, '  (n = ', num2str(currSummary{iCond, 1}), ')'];
    end
   
    % Plot figures
    dffCurrConds = combinedDffAvg{eventInd}(:,:,:, currConds);
    range = calc_range(dffCurrConds, rangeScalar, rangeType);
    [~, ~] = plot_heatmaps(dffCurrConds, analysisMetadata, range, plotTitles(currConds), sigma, 'fileName', fileName, 'makeVid', makeVid, 'saveDir', saveDir);
     
clear odorEventName eventInd currSummary currCondNames currConds sigma rangeType rangeScalar makeVid saveDir fileName plotTitles dffCurrConds range
catch(foldME); rethrow(foldME); end  
%% =================================================================================================
%   BEHAVIOR HEATMAPS                                   
%%%=================================================================================================
    %% PLOT OVERALL MEAN dF/F ACROSS BEHAVIORAL STATES
try
    behaviorNames = {'Locomotion', 'Grooming', 'IsoMove', 'AnyMove'};
    actionNum = [1];

    smoothingSigma = [0.6]; 
    rangeType = 'max';
    rangeScalar = 0.8;
    makeVid = 0;
    saveDir = [];
    fileName = ['All_frame_dFF_', behaviorNames{actionNum}, '_Heatmaps'];
    titleStr = {['dF/F - ', behaviorNames{actionNum}, ' vs. Quiescence']};
    
    % Load dF/F data
    load(fullfile(parentDir, ['actionDff_', behaviorNames{actionNum}, '.mat'])) % --> 'actionDff'

    % Calculate absolute max dF/F value across all planes and action states
    range = calc_range(actionDff, rangeScalar, rangeType);

    % Plot figures
    [f, ~] = plot_heatmaps(actionDff, analysisMetadata, range, titleStr, smoothingSigma, 'fileName', fileName, 'makeVid', makeVid, ...
                           'saveDir', saveDir);
clear smoothingSigma rangeType rangeScalar makeVid saveDir fileName titleStr range   
catch(foldME); rethrow(foldME); end
    %% PLOT INTERACTION HEATMAPS FOR SOME TRIAL CONDITIONS
try
    % Show summary again
    eventInd = ~cellfun(@isempty, strfind(primaryEventNames, 'locomotion'));
    currSummary = allCondSummaries{eventInd};
    disp(currSummary)
                 
    currCondNames = repmat(allCondNames{eventInd}, 2, 1);
    
    % Calculate absolute max dF/F value across all planes and stim types
    currConds = [1 2 3];
    sigma = [0.6]; 
    rangeType = 'Max';
    rangeScalar = 0.5;
    makeVid = 1;
    saveDir = [];
    fileName = 'Locomotion_Onset_Heatmaps';

    plotTitles = [];
    for iCond = currConds
        plotTitles{iCond} = [regexprep(currCondNames{iCond}, '_', '\\_'), ' ', currSummary{iCond, 3}{:}, '  (n = ', num2str(currSummary{iCond, 1}), ')'];
    end
   
    % Plot figures
    dffCurrConds = combinedDffAvg{eventInd}(:,:,:, currConds);
    range = calc_range(dffCurrConds, rangeScalar, rangeType);
    [~, ~] = plot_heatmaps(dffCurrConds, analysisMetadata, range, plotTitles(currConds), sigma, 'fileName', fileName, 'makeVid', makeVid, 'saveDir', saveDir); 

    clear eventInd currSummary currCondNames currConds sigma rangeType rangeScalar makeVid saveDir fileName plotTitles range

catch(foldME); rethrow(foldME); end

%% =================================================================================================
%           ROI-BASED ANALYSES                                   
%%%=================================================================================================
    %% PLOT AND SAVE NEW ROIs
try
ROIselectionGui();
catch(foldME); rethrow(foldME); end
    %% LOAD ROI DATA AND PLOT ROIs
try
parentDir = ['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid)];
metaDataFileName = 'ROI_metadata.mat';
dffDataFileName = 'ROI_Data_Avg.mat';

% Load metadata
load(fullfile(parentDir, metaDataFileName)); % --> ROImetadata(.mask, .xi, .yi, .plane, .color, .refImg)
analysisMetadata.ROImetadata = ROImetadata;
analysisMetadata.nROIs = numel(ROImetadata); nROIs = analysisMetadata.nROIs;

% Load imaging and dF/F data
load(fullfile(parentDir, dffDataFileName)); % --> 'ROIDataAvg', 'ROIDffAvg' ([volume, trial, ROI])
disp('ROI data loaded')

% Plot ROIs 
plot_ROIs(ROImetadata);

clear metaDataFileName dffDataFileName 
catch(foldME); rethrow(foldME); end
    %% PLOT 2D HEATMAPS OF ROI FLUORESCENCE THROUGHOUT EXPERIMENT
try
s = stimSepTrials;
saveDir = uigetdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');

useDff = 0;
% % 
trialGroups = [];
plotTitleSuffix = '';
fileNameSuffix = '_AllTrials';

% trialGroups = [s.OdorA + 2 * s.OdorB]; 
% plotTitleSuffix = ['  —  ', stimNames{1}, ' (top) vs. ', stimNames{2}, ' (bottom)']; %
% fileNameSuffix = ['_OdorAvsOdorB'];

%---------------------------------------------------------------------------------------------------

for iROI = 1:nROIs
    
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

%     if ~isempty(trialGroups)
%         trialGroups(:, ~goodTrials) = [];
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
    plot_2D_summary(analysisMetadata, flData(:,:,iROI)', ...
        'plotAxes', ax, ...
        'trialGroups', trialGroups, ...
        'titleStr', [regexprep(expDate, '_', '\\_'), '  -  ROI ', num2str(iROI), '   Raw fluorescence (AU)', plotTitleSuffix], ...
        'saveDir', saveDir, ...
        'colormap', cm, ...
        'fileName', fileName);
end%iROI
catch(foldME); rethrow(foldME); end
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
        [baselineData, respData] = extract_event_volumes(eventList, combFilterVecs{iType}(goodEvents,iCond), baselineDur, respDur, analysisMetadata, ...
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
catch(foldME); rethrow(foldME); end
            %% PLOT EVENT-ALIGNED dF/F WITHIN ROIs
try
% Show summary again
eventName = 'locomotion';
fileNamePrefix = 'Locomotion_offset_responses_';
shadeDur = 0;
eventInd = contains(primaryEventNames, eventName);
currSummary = allCondSummaries{eventInd};
disp(currSummary)

currConds = [4 5 6];
currCondNames = allCondNames{eventInd}(currConds);

currDffData = ROIEventDff{eventInd}(currConds); % --> {cond}[volume, event, ROI]

saveDir = 0;
saveDir = uigetdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');


ROIlist = 1:size(currDffData{1}, 3);
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
    imshow(analysisMetadata.refImg{currPlane}, [0 MAX_INTENSITY]);
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
catch(foldME); rethrow(foldME); end
    %% PLOT MEAN ROI dF/F THROUGHOUT TRIAL FOR ONE OR MORE STIM TYPES
try
    saveDir = uigetdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
    plotTitle = regexprep([expDate, ' - ', stimNames{1}, ' (top) ', stimNames{2}, ' (bottom)'], '(?<!\\)_', '\\_');
    fileNamePrefix = 'Whole_Trial_Responses_';
    s = analysisMetadata.stimSepTrials;
    eventShading = [10 13];
    filterVecs = logical([s.OdorA; s.OdorB] .* repmat(goodTrials, 2, 1));
    
    singleTrials = 1;
    singleTrialAlpha = 0.5;
    stdDevShading = 1;
    
    ROIlist = 1:size(ROIDffAvg, 3);
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
        imshow(analysisMetadata.refImg{currPlane}, [0 analysisMetadata.MAX_INTENSITY]);
        plot(xData, yData, 'Color', 'g');

        clear ax
        for iPlot = 1:nPlots
            currDffAvg = ROIDffAvg(:, filterVecs(iPlot, :), iROI); % --> [volume, trial]
            subaxis(nRows, 2, 1, iPlot+1, 2,1)  %( (iPlot*3 + 1):(iPlot*3 + 3) ))                                                     
            ax(iPlot) = gca;
            plot_ROI_data(ax(iPlot), currDffAvg, 'EventShading', eventShading, ...
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
catch(foldME); rethrow(foldME); end    
    %% PLOT MEAN ROI dF/F THROUGHOUT TRIAL FOR EARLY VS LATE TRIALS
try
    saveDir = uigetdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
    
    plotTitle = regexprep([expDate, ' - All stims - early vs late trials'], '_', '\\_');
    fileNamePrefix = 'EarlyVsLateTrials_';
    s = analysisMetadata.stimSepTrials.OdorA + analysisMetadata.stimSepTrials.OdorB;
    
    eventShading = [10 13];
    
    singleTrials = 1;
    singleTrialAlpha = 0.35;
    stdDevShading = 1;
    outlierSD = 5;
    legendStr = {'Trials 1:35', 'Trials 36:70', ''};
    
    trialGroups = ones(1, nTrials);
    trialGroups(36:end) = 2;
    trialGroups(91:end) = 3;
%     trialGroups(80:end) = 3;
    trialGroups = trialGroups(logical(s .* goodTrials));

    ROIlist = 1:size(ROIDffAvg, 3);
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
        imshow(analysisMetadata.refImg{currPlane}, [0 analysisMetadata.MAX_INTENSITY]);
        plot(xData, yData, 'Color', 'g');
%         title(['Plane #', num2str(analysisMetadata.ROImetadata{iROI}(1).plane)])

        % Create plot
        currDffAvg = ROIDffAvg(:,logical(s .* goodTrials), iROI); % --> [volume, trial]
        subaxis(2, 3, [4:6])
        ax = gca;
        plot_ROI_data(ax, currDffAvg, 'EventShading', eventShading, ...
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
catch(foldME); rethrow(foldME); end
    %% PLOT MEAN ROI dF/F THROUGHOUT TRIAL COLOR CODED BY BEHAVIOR
try
    saveDir = uigetdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');
    plotTitle = regexprep([expDate, ' - ', stimNames{1}, ' (top) ', stimNames{2}, ' (bottom)'], '(?<!\\)_', '\\_');
    fileNamePrefix = 'Behavior_Coded_Whole_Trial_Responses_';
    s = analysisMetadata.stimSepTrials;
    eventShading = [10 13];
    annotValues = [2 0];    
    
    trialFilterVecs = logical([s.OdorA; s.OdorB] .* repmat(goodTrials, 2, 1));
    
    singleTrials = 1;
    singleTrialAlpha = 0.4;
    stdDevShading = 0;
    
    ROIlist = 1:size(ROIDffAvg, 3);
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
        imshow(analysisMetadata.refImg{currPlane}, [0 analysisMetadata.MAX_INTENSITY]);
        plot(xData, yData, 'Color', 'g');
%         title(['Plane #', num2str(analysisMetadata.ROImetadata{iROI}(1).plane)])
        
        clear ax
        yL = [];
        for iPlot = 1:nPlots
            currDffAvg = ROIDffAvg(:, trialFilterVecs(iPlot, :), iROI); % --> [volume, trial]
            annotData = annotationTypes{6}.volAnnotArr;
            currAnnotArr = annotData(trialFilterVecs(iPlot, :), :);

            subaxis(nRows, 3, ( (iPlot*3 + 1):(iPlot*3 + 3) ))
            ax(iPlot) = gca; 
                        
            plot_ROI_data(ax(iPlot), currDffAvg, 'AnnotArray', currAnnotArr', ... 
                                                 'AnnotValues', annotValues', ...
                                                 'EventShading', eventShading, ...
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
catch(foldME); rethrow(foldME); end    
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
catch(foldME); rethrow(foldME); end

%% %=================================================================================================
%%% FICTRAC DATA PLOTTING
%% %================================================================================================
s = stimSepTrials;

currStim = s.OdorA + s.OdorB;
smWin = 3;
nBins = 50;
thresh = 0.5;
currTrial = 52;

try
    
% Pull out good trials from data 
mmSpeedData = ftData.moveSpeed(:,logical(goodTrials .* currStim)) * FRAME_RATE * 4.5;   % --> [frame, trial] (mm/sec)
dHD = rad2deg(ftData.yawSpeed(:,logical(goodTrials .* currStim)) * FRAME_RATE);         % --> [frame, trial] (deg/sec)
dHDSmooth = movmean(movmean(abs(dHD), smWin, 1), smWin, 1);                                 
% goodFl = ROIDffAvg(:,logical(goodTrials .* currStim),:);
goodFl = ROIDataAvg(:,logical(goodTrials .* currStim),:);   % --> [volume, trial, ROI]

% Normalize raw fluorescence for each ROI
for iROI = 1:size(goodFl, 3)
    currFl = goodFl(:,:,iROI);
    goodFl(:,:,iROI) = currFl ./ max(currFl(:));
end

goodTrialNums = 1:nTrials;
goodTrialNums(~logical(goodTrials .* currStim)) = [];
goodTrialTypes = analysisMetadata.trialType;
goodTrialTypes(~logical(goodTrials .* currStim)) = [];

% Scale fluorescence data so min is zero
goodFl = goodFl - min(goodFl(:));

% Calculate correlations between the fly's speed, yaw, and dF/F
allSpeed = as_vector(movmean(movmean(mmSpeedData(volFrames,:), smWin, 2), smWin, 2));
allYaw =  as_vector(dHDSmooth(volFrames,:)); % Note that this is currently directionless yaw speed

nROIs = 1%size(ROIDataAvg, 3);
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
    imshow(refImg{ROImetadata{iROI}(1).plane}, [0 analysisMetadata.MAX_INTENSITY]);hold on
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
if ~isempty(analysisMetadata.stimOnsetTimes)
    shadeTimes = [analysisMetadata.stimOnsetTimes(1), analysisMetadata.stimOnsetTimes(1) + analysisMetadata.stimDurs(1)];
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
if ~isempty(analysisMetadata.stimOnsetTimes)
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
if ~isempty(analysisMetadata.stimOnsetTimes)
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

catch(foldME); rethrow(foldME); end
%% =================================================================================================
%  Create vid with trial-by-trial FicTrac plots
%% =================================================================================================
try
saveDir = uigetdir(['B:\Dropbox (HMS)\2P Data\Imaging Data\', expDate, '\sid_', num2str(sid), '\Analysis'], 'Select a save directory');

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
    title([analysisMetadata.trialType{iTrial}, ' (', num2str(analysisMetadata.stimOnsetTimes(iTrial)), ...
          '-', num2str(analysisMetadata.stimOnsetTimes(iTrial) + analysisMetadata.stimDurs(iTrial)), ' sec)'])
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

catch(foldME); rethrow(foldME); end
%% =================================================================================================
%           OTHER ANALYSES                                   
%%%=================================================================================================

    %% PCA
try
% Pull out data for one plane
planeNum = 1;
planeData = pcaData(:,:,:,planeNum);

figure(2); clf;
subplot(2,2,1)
imshow(analysisMetadata.refImg{planeNum},[0 analysisMetadata.MAX_INTENSITY])
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
catch(foldME); rethrow(foldME); end
    %% VIEW RAW DATA FOR A SINGLE TRIAL AND PLANE
try
if ~exist('m', 'var')
    [dataFile, sessionDataPath, ~] = uigetfile('*.mat', 'Select an imaging data session file', 'B:\Dropbox (HMS)\2P Data\Imaging Data\');
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
catch(foldME); rethrow(foldME); end
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
catch(foldME); rethrow(foldME); end
    %% CREATE FICTRAC + BEHAVIOR VIDS
try
vidDir = 'B:\Dropbox (HMS)\2P Data\Behavior Vids\2018_06_29\_Movies';
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
goodTrialTypes = analysisMetadata.trialType;
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
    if ~isempty(analysisMetadata.stimOnsetTimes)
        onsetTimes = analysisMetadata.stimOnsetTimes;
        onsetTimes(~goodTrials) = [];
        stimDurs = analysisMetadata.stimDurs;
        stimDurs(~goodTrials) = [];
        shadeTimes = [analysisMetadata.stimOnsetTimes(iTrial), analysisMetadata.stimOnsetTimes(iTrial) + analysisMetadata.stimDurs(iTrial)];
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
        if ~isempty(analysisMetadata.stimOnsetTimes)
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
        if ~isempty(analysisMetadata.stimOnsetTimes)
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
catch(foldME); rethrow(foldME); end