
%% LOAD DATA

expDate = '2019_01_23_exp_1';
% expDate = '2018_12_18_exp_3';
sid = 0;

parentDir = fullfile('D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate);
savePath = fullfile('D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, ['sid_', num2str(sid)]);

annotFileName = '_Movies\autoAnnotations.mat';

try
% ----------------  Load stim metadata -------------------------------------------------------------
infoStruct = [];
stimDataFiles = dir(fullfile(parentDir, ['metadata*sid_', num2str(sid), '*.mat']));
infoStruct.stimOnsetTimes = []; infoStruct.stimDurs = []; infoStruct.trialType = []; infoStruct.outputData = [];
for iFile = 1:numel(stimDataFiles)
    
    % Load file    
    load(fullfile(parentDir, stimDataFiles(iFile).name)); % variable "metaData" with fields 'trialDuration', 'nTrials', 'stimTypes', 'sid', 'taskFile', 'outputData'
    
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
annotData = load(fullfile(parentDir, annotFileName)); % variables 'trialAnnotations', 'annotParams', 'ftData', 'flowArr', 'goodTrials', 'behaviorLabels', 'frameInfo'
annotData.nFrames = annotData.frameInfo.nFrames;
annotData.frameTimes = annotData.frameInfo.frameTimes;
annotData.FRAME_RATE = annotData.frameInfo.FRAME_RATE;
infoStruct = setstructfields(infoStruct, annotData);


% ----------------  Load FicTrac data --------------------------------------------------------------
ftData = load_fictrac_data(infoStruct, 'sid', sid, 'ParentDir', fullfile(parentDir, '\_Movies\FicTracData'));
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
% if ~isdir(fullfile(parentDir, '\_Movies\Analysis'))
%     mkdir(fullfile(parentDir, '\_Movies\Analysis'));
% end
if ~isdir(savePath)
    mkdir(savePath)
end
catch foldME; rethrow(foldME); end

%% SET UP PLOTTING VARIABLES

stimNames = {'EtOH\_neat'}; % {'EtOH\_neat', 'CO2\_2e-2'}; % 'Benzaldehyde\_e-1'
stimTrialGroups = s.OdorA; % [s.OdorA + 2 * s.OdorB]; %
stimGroupNames = {'OdorA'}; % {'OdorA', 'OdorB'}; % 
stimShading =  {[6 7; 8 11; 10 11]};% {[8 11]};%{[7 10 ; 10 13]};% {[4 6;10 13]};%
stimShadingColors = {'red', 'green', 'blue'}; % {'red', 'green'}; % 
rgbStimShadeColors = [rgb(stimShadingColors{1}); rgb(stimShadingColors{2}); rgb(stimShadingColors{3})];
groupBounds = [1:40:nTrials]; groupBounds(2:end) = groupBounds(2:end) - 1;
% groupBounds = [1, 40:60:nTrials-1]; groupBounds(2:end) = groupBounds(2:end) -1;
% groupBounds = [1, 60];
blockNames = {'OdorOnly', 'FW', 'OdorOnly2', 'BW', 'OdorOnly3', '20%'};

%% 2D BEHAVIOR SUMMARY
saveDir = uigetdir(savePath, 'Select a save directory');
% 
% trialGroups = [];
% plotTitleSuffix = '';
% fileNameSuffix = '_AllTrials';
% % '
% % % 
% trialGroups = stimTrialGroups; 
% plotTitleSuffix = make_plotTitleSuffix(stimNames); %
% fileNameSuffix = make_fileNameSuffix(stimGroupNames);

% % GROUP BY BLOCK TYPE
% groupNames = {'Odor only', 'Odor + photostim'};
% % groupNames = {'Odor + photostim', 'Odor only'};
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = ~mod(iBound, 2);
% end
% trialGroups(groupBounds(end):end) = ~mod(iBound + 1, 2);
% trialGroups = trialGroups + 1;
% plotTitleSuffix = make_plotTitleSuffix(groupNames);
% fileNameSuffix = '_PhotostimVsOdorOnly';
% 
% % GROUP CHRONOLOGICALLY BY BLOCKS
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(groupBounds(end):end) = iBound + 1;
% plotTitleSuffix = '';
% fileNameSuffix = '_Blocks_Separated';

try
    
% Create plot titles
plotName = 'Behavior Annotation';
titleString = [regexprep(expDate, '_', '\\_'), '    ', [plotName, ' summary ', plotTitleSuffix]]; % regex to add escape characters
annotArr = behaviorAnnotArr;

% Create figure
f = figure(7);clf
f.Position = [-1050 45 1020 950];
f.Color = [1 1 1];
ax = gca();
[~, ax, ~] = plot_behavior_summary_2D(infoStruct, annotArr, ax, titleString, trialGroups);
ax.FontSize = 14;
ax.Title.FontSize = 12;
ax.XLabel.FontSize = 14;

% Plot stim times
hold on
[nStimEpochs, idx] = max(cellfun(@size, stimShading, repmat({1}, 1, numel(stimShading))));
for iStim = 1:nStimEpochs
    stimOnsetFrame = stimShading{idx}(iStim, 1) * FRAME_RATE;
    stimOffsetFrame = stimShading{idx}(iStim, 2) * FRAME_RATE;
    plot(ax, [stimOnsetFrame, stimOnsetFrame], ylim(), 'Color', 'k', 'LineWidth', 2)
    plot(ax, [stimOffsetFrame, stimOffsetFrame], ylim(), 'Color', 'k', 'LineWidth', 2)
end

if saveDir
    fileName = ['2D_Annotation_Summary', fileNameSuffix, '_', ...
                regexprep(expDate, {'_', 'exp'}, {'', '_'})];
    export_fig(fullfile(saveDir, fileName), '-png', f);
    if ~isdir(fullfile(saveDir, 'figFiles'))
        mkdir(fullfile(saveDir, 'figFiles'))
    end
    savefig(f, fullfile(saveDir, 'figFiles', fileName));
end
catch foldME; rethrow(foldME); end

%% 1D BEHAVIOR SUMMARY
actionNames = {'NA', 'IsoMovement', 'Locomotion'};
saveDir = uigetdir(savePath, 'Select a save directory');
actionNum = [3]; % locomotionLabel = 3; noActionLabel = 0; isoMovementLabel = 1;
actionName = actionNames{actionNum};
figTitle = regexprep([expDate, '  —  Fly ', actionName, ' throughout trial'], '_', '\\_');
cm = [];
% 
% % ALL TRIALS
% trialGroups = [goodTrials];
% fileNameSuffix = ['_AllTrials_', actionName];
% groupNames = {'All trials'};

% % GROUP BY STIM TYPE
% trialGroups = stimTrialGroups .* goodTrials;
% fileNameSuffix = [make_fileNameSuffix(stimGroupNames), '_', actionName]; 
% groupNames = stimNames;

% % GROUP BY SINGLE BLOCKS
% groupNames = [];
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
%    groupNames{iBound} = ['Trials ', num2str(groupBounds(iBound)), '-', num2str(groupBounds(iBound + 1))];
% end
% trialGroups(groupBounds(end):end) = numel(groupBounds);
% groupNames = {'Odor only', 'Odor + photostim'};
% fileNameSuffix = '_SingleBlocks';
% cm = repmat([rgb('blue'); rgb('red')], 8, 1);
% % 
% % GROUP BY BLOCK TYPE
% groupNames = {'Odor only', 'Odor + photostim'};
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = ~mod(iBound, 2);
% end
% trialGroups(groupBounds(end):end) = ~mod(iBound + 1, 2);
% trialGroups = trialGroups + 1;
% fileNameSuffix = '_PhotostimVsOdorOnly';
% % 
% % PLOT AND COLOR EACH BLOCK SEPARATELY
% groupNames = blockNames;
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(groupBounds(end):end) = numel(groupBounds);
% fileNameSuffix = '_Blocks_Separated';



try
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
        cm = [rgb('blue'); rgb('red'); rgb('green'); rgb('magenta'); rgb('cyan'); rgb('gold'); rgb('lime')];
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
    
    legend(groupNames, 'FontSize', 14, 'Location', 'Best', 'AutoUpdate', 'off')
    ax.XLim = [20 nFrames-20]; % to improve plot appearance
    ax.YLim = [0 1];
    suptitle(figTitle);
end

% Add shading during stimulus presentations
[nStimEpochs, idx] = max(cellfun(@size, stimShading, repmat({1}, 1, numel(stimShading))));
yL = ylim();
for iStim = 1:size(stimShading{idx}, 1)
    stimOnset = stimShading{idx}(iStim, 1) * FRAME_RATE;
    stimOffset = stimShading{idx}(iStim, 2) * FRAME_RATE;
    plot_stim_shading([stimOnset, stimOffset], 'Color', rgb(stimShadingColors{iStim}))
end

if saveDir
    fileName = ['1D_Annotation_Summary', fileNameSuffix, '_', ...
                regexprep(expDate, {'_', 'exp'}, {'', '_'})];
    export_fig(fullfile(saveDir, fileName), '-png', f);
    if ~isdir(fullfile(saveDir, 'figFiles'))
        mkdir(fullfile(saveDir, 'figFiles'))
    end
    savefig(f, fullfile(saveDir, 'figFiles', fileName));
end
catch foldME; rethrow(foldME); end

%% 2D FICTRAC SUMMARY

saveDir = uigetdir(savePath, 'Select a save directory');
fontSize = 12;

ftVarName = 'moveSpeed'; % 'moveSpeed', 'fwSpeed', 'yawSpeed'%
sdCap = 3.5;
smWin = 9;
cmName = @parula;
figTitle = [regexprep(expDate, '_', '\\_'), '  —  FicTrac ', ftVarName];


% % ALL TRIALS
% trialGroups = [];
% fileNameSuffix = ['_AllTrials'];
% plotTitleSuffix = '';

% % % % % % % % % 
% % GROUP BY STIM TYPE
% trialGroups = stimTrialGroups .* goodTrials; 
% fileNameSuffix = make_fileNameSuffix(stimGroupNames);
% plotTitleSuffix = make_plotTitleSuffix(stimNames);
% % 
% % GROUP BY BLOCK TYPE
% groupNames = {'Odor only', 'Odor + photostim'};
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = ~mod(iBound, 2);
% end
% trialGroups(groupBounds(end):end) = ~mod(iBound + 1, 2);
% trialGroups = (trialGroups + 1) .* goodTrials;
% plotTitleSuffix = make_plotTitleSuffix(groupNames);
% fileNameSuffix = '_PhotostimVsOdorOnly';
% 
% % GROUP BY BLOCKS
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%     trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(groupBounds(end):end) = iBound + 1;
% trialGroups = trialGroups .* goodTrials;
% plotTitleSuffix = '';
% fileNameSuffix = '_Blocks_Separated';


try

% Extract FicTrac data
rawData = infoStruct.ftData.(ftVarName);          % --> [frame, trial]
rawData = rawData';                                     % --> [trial, frame]
if strcmp(ftVarName, 'yawSpeed')
    plotData = abs(rad2deg(rawData .* FRAME_RATE));    	% --> [trial, frame] (deg/sec)
    figTitle = [figTitle, ' (deg/sec)'];
else
    plotData = rawData .* FRAME_RATE .* 4.5;            % --> [trial, frame] (mm/sec)
    figTitle = [figTitle, ' (mm/sec)'];
end

% Cap values at n SD above mean
capVal = mean(plotData(:), 'omitnan') + (sdCap * std(plotData(:), 'omitnan'));
plotData(plotData > capVal) = capVal;

% Smooth data
smPlotData = movmean(plotData, smWin, 2);

% Create colormap
cm = cmName(numel(unique(smPlotData)));
if ~isempty(trialGroups)
    cm = [0 0 0; cm(2:end, :)];
end

% Plot data
titleStr = [figTitle, plotTitleSuffix];
[~, ax, f] = plot_2D_summary(infoStruct, smPlotData, ...
                'trialGroups', trialGroups, ...
                'titleStr', titleStr, ...
                'sampRate', FRAME_RATE, ...
                'colormap', cm ...
                );
ax.Title.FontSize = fontSize;

% Plot stim times
colorbar
hold on

% Plot stim times
hold on
[nStimEpochs, idx] = max(cellfun(@size, stimShading, repmat({1}, 1, numel(stimShading))));
for iStim = 1:nStimEpochs
    stimOnsetFrame = stimShading{idx}(iStim, 1) * FRAME_RATE;
    stimOffsetFrame = stimShading{idx}(iStim, 2) * FRAME_RATE;
    plot(ax, [stimOnsetFrame, stimOnsetFrame], ylim(), 'Color', 'k', 'LineWidth', 2)
    plot(ax, [stimOffsetFrame, stimOffsetFrame], ylim(), 'Color', 'k', 'LineWidth', 2)
end

if saveDir
    fileName = ['2D_FicTrac_', ftVarName, '_Summary', fileNameSuffix, '_', ...
                regexprep(expDate, {'_', 'exp'}, {'', '_'})];
    export_fig(fullfile(saveDir, fileName), '-png', f);
    if ~isdir(fullfile(saveDir, 'figFiles'))
        mkdir(fullfile(saveDir, 'figFiles'))
    end
    savefig(f, fullfile(saveDir, 'figFiles', fileName));
end

catch foldME; rethrow(foldME); end

%% 1D TRIAL-AVGERAGED FICTRAC SUMMARY

saveDir = uigetdir(savePath, 'Select a save directory');

includeQuiescence = 0;
if ~includeQuiescence
    fileNameSuffix = 'NoQuiescence_';
else
    fileNameSuffix = '';
end
figTitle = [expDate, '  —  Trial-Averaged FicTrac data'];
trialGroups = goodTrials';
smWin = 11;
cm = [];


% % ALL TRIALS
% trialGroups = [goodTrials];
% fileNameSuffix = [fileNameSuffix, 'AllTrials'];
% groupNames = {'All trials'};


% % GROUP BY STIM TYPE
% trialGroups = stimTrialGroups .* goodTrials; 
% fileNameSuffix = [fileNameSuffix, 'StimTypeComparison']; 
% groupNames = stimNames;

% % % 
% % GROUP BY SINGLE BLOCKS
% groupNames = [];
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
%    groupNames{iBound} = ['Trials ', num2str(groupBounds(iBound)), '-', num2str(groupBounds(iBound + 1))];
% end
% trialGroups(groupBounds(end):end) = numel(groupBounds);
% groupNames = {'Odor only', 'Odor + photostim'};
% fileNameSuffix = [fileNameSuffix, 'SingleBlocks'];
% cm = repmat([rgb('blue'); rgb('red')], 8, 1);
% %  
% % 
% % GROUP BY BLOCK TYPE
% groupNames = {'Odor only', 'Odor + photostim'};
% % groupNames = {'Odor only', 'Odor + photostim (10%)', 'Odor + photostim (25%)'};
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = ~mod(iBound, 2);
% end
% trialGroups(groupBounds(end):end) = ~mod(iBound + 1, 2);
% trialGroups = trialGroups + 1;
% fileNameSuffix = [fileNameSuffix, 'PhotoStimVsOdorOnly'];

% % PLOT AND COLOR EACH BLOCK SEPARATELY
% groupNames = blockNames;
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(groupBounds(end):end) = numel(groupBounds);
% fileNameSuffix = [fileNameSuffix, 'Blocks_Separated'];

try
    
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
axMoveSpeed = subaxis(3,1,1, 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.05); hold on
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
%     currFWSpeed = fwSpeed(:, trialGroups == iGroup);
    currYawSpeed = dHD(:, trialGroups == iGroup);
    currYawVel = yawVel(:, trialGroups == iGroup);
    
    if ~includeQuiescence
        currAnnotData = behaviorAnnotArr';
        currAnnotData(:, trialGroups ~= iGroup) = [];
        currXYSpeed(currAnnotData == 0) = nan;
%         currFWSpeed(currAnnotData == 0) = nan;
        currYawSpeed(currAnnotData == 0) = nan;
        currYawVel(currAnnotData == 0) = nan;
    end
    
    % Omit outliers
    outlierCalc = @(x) mean(x) + 4 * std(x);
    currXYSpeed(currXYSpeed >= outlierCalc(mmSpeedData(:))) = nan;
%     currFWSpeed(currFWSpeed >= outlierCalc(fwSpeed(:))) = nan;
    currYawSpeed(currYawSpeed >= outlierCalc(dHD(:))) = nan;
    currYawVel(currYawVel >= outlierCalc(yawVel(:))) = nan;

    meanSpeed = smooth(mean(currXYSpeed, 2, 'omitnan'), smWin);
%     meanFWSpeed = smooth(mean(currFWSpeed, 2, 'omitnan'), smWin);
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
axMoveSpeed.XTick = xTickFR;
axMoveSpeed.XTickLabel = xTickLabels;
axMoveSpeed.YLabel.String = 'XY Speed (mm/sec)';
axMoveSpeed.FontSize = 14;
legend(axMoveSpeed, groupNames, 'FontSize', 12, 'Location', 'NW', 'AutoUpdate', 'off')
axMoveSpeed.XLim = [9 nFrames-5]; % to improve plot appearance
if ~isempty(stimShading)
    [nStimEpochs, idx] = max(cellfun(@size, stimShading, repmat({1}, 1, numel(stimShading))));
    for iStim = 1:size(stimShading{idx}, 1)
        stimOnset = stimShading{idx}(iStim, 1) * FRAME_RATE;
        stimOffset = stimShading{idx}(iStim, 2) * FRAME_RATE;
        plot_stim_shading([stimOnset, stimOffset], 'Color', rgb(stimShadingColors{iStim}), 'Axes', ...
            axMoveSpeed);
    end
end

axYawVel.XTick = xTickFR;
axYawVel.XTickLabel = xTickLabels;
axYawVel.YLabel.String = 'Yaw Vel (CW = +)';
axYawVel.FontSize = 14;
legend(axYawVel, groupNames, 'FontSize', 12, 'Location', 'NW', 'AutoUpdate', 'off')
axYawVel.XLim = [9 nFrames-5]; % to improve plot appearance
if ~isempty(stimShading)
    [nStimEpochs, idx] = max(cellfun(@size, stimShading, repmat({1}, 1, numel(stimShading))));
    for iStim = 1:size(stimShading{idx}, 1)
        stimOnset = stimShading{idx}(iStim, 1) * FRAME_RATE;
        stimOffset = stimShading{idx}(iStim, 2) * FRAME_RATE;
        plot_stim_shading([stimOnset, stimOffset], 'Color', rgb(stimShadingColors{iStim}), 'Axes', ...
            axYawVel);
    end
end

axYawSpeed.XTick = xTickFR;
axYawSpeed.XTickLabel = xTickLabels;
axYawSpeed.YLabel.String = 'Yaw Speed (deg/sec)';
axYawSpeed.FontSize = 14;
legend(axYawSpeed, groupNames, 'FontSize', 12, 'Location', 'NW', 'AutoUpdate', 'off')
axYawSpeed.XLim = [9 nFrames-5]; % to improve plot appearance
if ~isempty(stimShading)
    [nStimEpochs, idx] = max(cellfun(@size, stimShading, repmat({1}, 1, numel(stimShading))));
    for iStim = 1:size(stimShading{idx}, 1)
        stimOnset = stimShading{idx}(iStim, 1) * FRAME_RATE;
        stimOffset = stimShading{idx}(iStim, 2) * FRAME_RATE;
        plot_stim_shading([stimOnset, stimOffset], 'Color', rgb(stimShadingColors{iStim}), 'Axes', ...
            axYawSpeed);
    end
end

suptitle(regexprep([figTitle, '  —  ', fileNameSuffix], '_', '\\_'));
if saveDir

    % Create filename
    fileName = regexprep(['1D_FicTrac_Summary_', fileNameSuffix, '_', ...
                            regexprep(expDate, {'_', 'exp'}, {'', '_'})], '_', '\_');
    
    % Save figure files
        export_fig(fullfile(saveDir, fileName), '-png', f);
        if ~isdir(fullfile(saveDir, 'figFiles'))
            mkdir(fullfile(saveDir, 'figFiles'))
        end
        savefig(f, fullfile(saveDir, 'figFiles', fileName));
end%if

catch foldME; rethrow(foldME); end


%% VOLUME-AVGERAGED FICTRAC SUMMARY

saveDir = uigetdir(savePath, 'Select a save directory');
smWin = 11;
sdCap = 3.5;
figTitle = [expDate, '  —  Volume-Averaged FicTrac data'];
cm = [];

% % ALL TRIALS
% trialGroups = [goodTrials];
% fileNameSuffix = [fileNameSuffix, 'AllTrials'];
% groupNames = {'All trials'};

% % GROUP BY STIM TYPE
% trialGroups = stimTrialGroups .* goodTrials; 
% fileNameSuffix = [fileNameSuffix, 'StimTypeComparison']; 
% groupNames = stimNames;

% PLOT AND COLOR EACH BLOCK SEPARATELY
groupNames = blockNames;
trialGroups = zeros(1, nTrials);
for iBound = 1:numel(groupBounds)-1
   trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
end
trialGroups(groupBounds(end):end) = numel(groupBounds);
fileNameSuffix = '_Blocks_Separated';

try

% Extract FicTrac data
mmSpeedData = ftData.moveSpeed * FRAME_RATE * 4.5;  % --> [frame, trial] (mm/sec)
dHD = abs(rad2deg(ftData.yawSpeed * FRAME_RATE));        % --> [frame, trial] (deg/sec)
fwSpeed = ftData.fwSpeed * FRAME_RATE * 4.5;        % --> [frame, trial  (mm/sec)
allFtData = cat(3, mmSpeedData, fwSpeed, dHD);
nFrames = size(mmSpeedData, 1);

% Cap values at n SD above mean
for iVar = 1:size(allFtData, 3)
    currData = allFtData(:,:,iVar);
    capVal = mean(currData(:), 'omitnan') + (sdCap * std(currData(:), 'omitnan'));
    currData(currData > capVal) = capVal;
    allFtData(:,:,iVar) = currData;
end

% Smooth data
allFtData = movmean(allFtData, 3, 2);

% smXYData = movmean(mmSpeedData, smWin, 2);
% smYawData = movmean(dHD, smWin, 2);
% smFWData = movmean(fwSpeed, smWin, 2);

% Create figure
f = figure(4); clf; hold on
f.Position = [-1100 50 900 930];
f.Color = [1 1 1];

% Create axes
M = 0.02;
P = 0.00;
axMoveSpeed = subaxis(3,1,1, 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.05); hold on
axFWSpeed = subaxis(3,1,2, 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.06); hold on
axYawSpeed = subaxis(3,1,3, 'S', 0, 'M', M, 'PB', 0.06, 'PL', 0.06); hold on

% Separate trialGroups
nGroups = length(unique(trialGroups(trialGroups ~= 0)));
if isempty(cm)
    cm = parula(nGroups);
    cm = [rgb('blue'); rgb('red'); rgb('green'); rgb('magenta'); rgb('cyan'); rgb('gold'); rgb('lime')];
end

groupArr = []; groupSize = [];
for iGroup = 1:nGroups
    if iGroup == 1
        groupArr = allFtData(:, trialGroups == iGroup, :);
        groupSize(1) = size(groupArr, 2);
    else
        groupArr = cat(2, groupArr, allFtData(:, trialGroups == iGroup, :));
        groupSize(end + 1) = size(groupArr, 2);
    end
end
plotArr = squeeze(mean(groupArr, 1)); % --> [trial, var]

% Plot data
axes(axMoveSpeed)
plot(plotArr(:,1), 'linewidth', 2)
for iGroup = 1:nGroups-1
    yL = ylim;
    plot([groupSize(iGroup), groupSize(iGroup)], [0 yL(2)], 'LineWidth', 2, 'Color', 'r')
end
xlim([0, size(plotArr, 1)])
% ylim([0 yL(2)])

axes(axFWSpeed)
plot(plotArr(:,2), 'linewidth', 2)
for iGroup = 1:nGroups-1
    yL = ylim;
    plot([groupSize(iGroup), groupSize(iGroup)], [0 yL(2)], 'LineWidth', 2, 'Color', 'r')
end
xlim([0, size(plotArr, 1)])
% ylim([0 yL(2)])

axes(axYawSpeed)
plot(plotArr(:,3), 'linewidth', 2)
for iGroup = 1:nGroups-1
    yL = ylim;
    plot([groupSize(iGroup), groupSize(iGroup)], [0 yL(2)], 'LineWidth', 2, 'Color', 'r')
end
xlim([0, size(plotArr, 1)])
% ylim([0 yL(2)])

catch foldME; rethrow(foldME); end


%% PLOT OVERLAID 2D MOVEMENT DATA

startTime = 14;
plotLen = 1;
limScalar = 0.5;
alpha = 0.4;
showMean = 1;
saveFig = 0;

s = stimSepTrials;
shadeFrames = {round(stimShading{:} * FRAME_RATE)};
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

% % % % % % 
% % GROUP BY BLOCK TYPE
% groupNames = {'Odor only', 'Odor + photostim'};
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = ~mod(iBound, 2);
% end
% trialGroups(groupBounds(end):end) = ~mod(iBound + 1, 2);
% trialGroups = trialGroups + 1;
% fileNameSuffix = [fileNameSuffix, 'PhotoStimVsOdorOnly'];

% % GROUP BY ALTERNATING BLOCKS
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(groupBounds(end):end) = iBound + 1;
% plotTitleSuffix = '';
% fileNameSuffix = '_Blocks_Separated';
% cm = jet(numel(unique(trialGroups)));

% % PLOT AND COLOR EACH BLOCK SEPARATELY
% groupNames = blockNames;
% trialGroups = zeros(1, nTrials);
% for iBound = 1:numel(groupBounds)-1
%    trialGroups(groupBounds(iBound):groupBounds(iBound + 1)) = iBound;
% end
% trialGroups(groupBounds(end):end) = numel(groupBounds);
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
    % Create filename
    fileName = regexprep(['Movement_Overlay_', num2str(startTime), '-', num2str(startTime + plotLen)...
        '_sec', fileNameSuffix, '_', ...
        regexprep(expDate, {'_', 'exp'}, {'', '_'})], '_', '\_');
    export_fig(fullfile(saveDir, fileName), '-png', f);
    if ~isdir(fullfile(saveDir, 'figFiles'))
        mkdir(fullfile(saveDir, 'figFiles'))
    end
    savefig(f, fullfile(saveDir, 'figFiles', fileName));
end

catch foldME; rethrow(foldME); end






