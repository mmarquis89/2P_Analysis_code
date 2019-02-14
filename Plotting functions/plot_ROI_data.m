function ax = plot_ROI_data(ax, ROIDffAvg, varargin)
%=========================================================================================================
% PLOT MEAN dF/F DATA WITHIN ROI THROUGHOUT TRIAL
%
% Creates a plot of trial-averaged dF/F data for an entire trial in a specified set of axes. There
% are various optional additions such as single trial plotting and standard deviation/stimulus event
% shading. Data is slightly smoothed by default but this can be disabled.
%
% INPUTS:
%
%       ax          = the handle to the axes where the plot should be created
%
%       ROIDffAvg   = the array of dF/F data with format: [volume, trial]
%
%
% OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%
%       'TrialGroups'      = (default: []) 1 x nTrials numeric array specifying two or more trial groups
%                            (numbers starting at one) to split the data into. Each group will be assigned
%                            a color (the individual trials will no longer be plotted chronologically),
%                            and instead of plotting the average of all trials each group will
%                            have its own trial-averaging line.
%
%       'OutlierSD'        = (default: 5) the number of standard deviations to use as an outlier
%                            exclusion threshold.
%
%       'SingleTrials'     = (default: 1) a boolean specifying whether to plot all individual trials
%                            in the background of the mean dF/F line.
%
%       'StdDevShading'    = (default: 1) boolean specifying whether to shade 1 standard deviation
%                            above and below the mean dF/F line.
%
%       'EventShading'     = (default: []) n x 2 element vector specifying times in seconds to shade to
%                            indicate stimulus presentation
%
%       'EventShadeColor'  = (default: [0 0 0]) the RGB color of each event shading 
%
%       'SmoothWinSize'    = (default: 3) the width in volumes of the moving average smoothing window
%                            that will be applied to BOTH the mean dF/F plot and the individual trial
%                            lines. Use a window size of "1" to skip smoothing entirely.
%
%       'VolOffset'        = (default: 0) a scalar number to add to the volumes for the purpose of
%                            labeling the X-axis. For example, if you are imaging at 10 volumes/sec
%                            and want your X-axis labeling to start at -2 seconds, use -20.
%
%       'SingleTrialAlpha' = (default: 1) a value from 0-1 specifying the transparency of the single
%                            trial lines.
%
%       'Legend'           = (default: []) a cell array of strings with one legend entry for each
%                            trialGroup you have provided
%
%       'VolumeRate'       = (default: 6.44) the rate (in hz) at which imaging volumes were acquired
%
%       'AnnotArray'       = (default: []) an array of behavior annotation data with the same dimensions
%                            as ROIDffAvg. If one is provided and SingleTrials == 1, the individual trials
%                            will be plotted in segments that are color-coded by behavior annotation in
%                            each volume. Note that the trials should not be divided into multiple groups
%                            if using this option or the coloring of the results will be confusing.
%
%       'AnnotValues'      = (default:[]) a vector of numbers to include when plotting behavior annotation
%                            segments. By default, all segments will be plotted, but passing a restricted
%                            list will only plot segments when their annotation value is a member of this 
%                            list.
%       'Colormap'         = (default: jet(numel(unique(annotArr))) )
%
%=========================================================================================================

% Parse optional arguments
p = inputParser;
addParameter(p, 'TrialGroups', []);
addParameter(p, 'OutlierSD', 5);
addParameter(p, 'SingleTrials', 1);
addParameter(p, 'StdDevShading', 1);
addParameter(p, 'EventShading', []);
addParameter(p, 'EventShadeColor', [0 0 0]);
addParameter(p, 'SmoothWinSize', 3);
addParameter(p, 'VolOffset', 0);
addParameter(p, 'SingleTrialAlpha', 1);
addParameter(p, 'Legend', []);
addParameter(p, 'VolumeRate', 6.44);
addParameter(p, 'AnnotArray', []);
addParameter(p, 'AnnotValues', []);
addParameter(p, 'Colormap', []);
parse(p, varargin{:});
trialGroups = p.Results.TrialGroups;
outlierSD = p.Results.OutlierSD;
singleTrials = p.Results.SingleTrials;
shadeSDs = p.Results.StdDevShading;
eventShading = p.Results.EventShading;
eventShadeColor = p.Results.EventShadeColor;
smoothWin = p.Results.SmoothWinSize;
volumeRate = p.Results.VolumeRate;
volOffset = p.Results.VolOffset;
singleTrialAlpha = p.Results.SingleTrialAlpha;
meanLineLegend = p.Results.Legend;
annotArr = p.Results.AnnotArray;
annotValues = p.Results.AnnotValues;
cm = p.Results.Colormap;

% Setup variables
nVolumes = size(ROIDffAvg, 1);
volTimes = ([1:1:nVolumes] + volOffset) ./ volumeRate;
nTrials = size(ROIDffAvg, 2);
if ~isempty(trialGroups)
    nGroups = numel(unique(trialGroups(trialGroups ~= 0)));
else
    nGroups = 1;
    trialGroups = ones(1, nTrials);
end
if isempty(annotValues)
    annotValues = unique(annotArr(:));
end

% Format axes
hold on
ax.YLabel.String =  'dF/F';
ax.XLabel.String = 'Time (s)';
xlim(ax, [min(volTimes), max(volTimes)]);

% % Discard any trials that are too many SDs from mean
% outliers = zeros(1, size(ROIDffAvg, 2));
% for iTrial = 1:size(ROIDffAvg, 2)
%     if max(abs(ROIDffAvg(:, iTrial) - trialAvgDff) > (outlierSD * mean(std(ROIDffAvg, 0, 2))))
%         outliers(iTrial) = 1;
%         
%     end
% end
% if sum(outliers) > 0
%     disp(['Omitting ' num2str(sum(outliers)), ' outlier trials'])
% end
% ROIDffAvg(:, logical(outliers)) = [];
% trialGroups(logical(outliers)) = [];
% if ~isempty(annotArr)
%     annotArr(:, logical(outliers)) = [];
% end

% Create colormap
if ~isempty(annotArr)
    if isempty(cm)
        cm =[rgb('RoyalBlue') * 0.6; ...
            rgb('Orange'); ...
            rgb('Crimson'); ...
            rgb('Cyan') * 0.85; ...
            rgb('Gold') * 0.9; ...
            rgb('Yellow'); ...
            rgb('Green'); ...
            rgb('Red'); ...
            rgb('blue'); ...
            rgb('Gold'); ...
            rgb('black') ...
            ];
    else
        cm = parula(numel(unique(annotArr)));
    end
elseif nGroups == 1
    cm = jet(nTrials);
else
    cm = [rgb('Blue'); ...
        rgb('Green'); ...
        rgb('Red'); ...
        rgb('Magenta'); ...
        rgb('Cyan'); ...
        rgb('Gold'); ...
        rgb('DarkRed'); ...
        rgb('Yellow'); ...
        rgb('Lime') ...
        ];
end

meanPlots = [];
for iGroup = 1:nGroups
    
    % Select plotting colors
    if nGroups == 1
        shadeColor = [0 0 0];
        meanLineColor = [0 0 0];
    else
        shadeColor = cm(iGroup, :);
        meanLineColor = shadeColor;
    end
    
    % Separate data from current trial group and calculate average dF/F
    groupAvgDff = mean(ROIDffAvg(:, trialGroups == iGroup), 2);
    groupDff = ROIDffAvg(:, trialGroups == iGroup); % --> [volume, trial]
    
    % Discard any trials that are too many SDs from mean
    outliers = zeros(1, size(ROIDffAvg, 2));
    for iTrial = 1:size(ROIDffAvg, 2)
        if trialGroups(iTrial) == iGroup
            if max(abs(ROIDffAvg(:, iTrial) - groupAvgDff) > (outlierSD * mean(std(ROIDffAvg, 0, 2))))
                outliers(iTrial) = 1;
            end
        end
    end
    groupOutliers = outliers(trialGroups == iGroup);
    if sum(groupOutliers) > 0
        disp(['Omitting ' num2str(sum(groupOutliers)), ' outlier trials'])
    end
    groupDff(:, logical(groupOutliers)) = [];
    if ~isempty(annotArr)
        annotArr(:, logical(groupOutliers)) = [];
    end
    groupStdDev = std(groupDff, 0, 2);
    
    % Re-calculate average without outliers
    groupAvgDff = mean(ROIDffAvg(:,  logical((trialGroups == iGroup) .* ~logical(outliers))), 2);

    % Plot individual trials in background
    if singleTrials
        if ~isempty(annotArr)
            uniqueAnnotVals = unique(annotArr(:));
            meanAnnotData = nan(size(groupDff, 1), size(groupDff, 2), length(uniqueAnnotVals)); % --> [volume, trial, annotVal]
        end
        for iTrial = 1:size(groupDff, 2)
            currData = smooth(groupDff(:, iTrial), smoothWin);
            
            if ~isempty(annotArr)
                
                % Plot trial dF/F in segments if coloring by annotation data
                currAnnotData = annotArr(:, iTrial);
                
                % Find index of each behavior transition in the current trial
                currAnnotDataStr = regexprep(num2str(currAnnotData'), ' ', '');
                postStartTransVols = regexp(currAnnotDataStr, '.(?=(0[234])|2[034]|3[024]|4[023])') + 1;
                startTransVols = regexp(currAnnotDataStr, '^0[234]'); % Necessary if there is a transition between volumes #1 and #2
                endTransVols = regexp(currAnnotDataStr, '[234]0$'); % Same but for the last few volumes
                transVols = [startTransVols, postStartTransVols, endTransVols];
                volTypes = currAnnotDataStr(transVols);
                
                % Plot each line segment in the proper color
                for iSeg = 1:numel(transVols)
                    if iSeg == 1
                        startVol = 1;
                        endVol = transVols(1);
                    else
                        startVol = transVols(iSeg - 1) + 1;
                        endVol = transVols(iSeg);
                    end
                    currVolType = str2double(volTypes(iSeg));
                    currColor = cm(currVolType + 1, :);
                    if ismember(currVolType, annotValues)
                        plt = plot(ax, volTimes(startVol:endVol+1), currData(startVol:endVol+1), 'color', currColor, 'LineWidth', 1);
                        plt.Color(4) = singleTrialAlpha;
                        meanAnnotData(startVol:endVol+1, iTrial, uniqueAnnotVals == currVolType) = currData(startVol:endVol+1); % --> [volume, trial, annotVal]
                    end
                    
                end
                
                if ~isempty(transVols)
                    % Still have to plot the final segment
                    startVol = transVols(end) + 1;
                    endVol = numel(currAnnotDataStr);
                    currVolType = str2double(currAnnotDataStr(end));
                    currColor = cm(currVolType + 1, :);
                    if ismember(currVolType, annotValues)
                        plt = plot(ax, volTimes(startVol:endVol), currData(startVol:endVol), 'color', currColor, 'LineWidth', 1);
                        plt.Color(4) = singleTrialAlpha;
                        meanAnnotData(startVol:endVol, iTrial, uniqueAnnotVals == currVolType) = currData(startVol:endVol); % --> [volume, trial, annotVal]
                    end
                else
                    % In case there's no behavior in the trial, just plot a single line
                    startVol = 1;
                    endVol = numel(currAnnotDataStr);
                    currVolType = str2double(currAnnotDataStr(1));
                    currColor = cm(currVolType + 1, :);
                    if ismember(currVolType, annotValues)
                        plt = plot(ax, volTimes(startVol:endVol), currData(startVol:endVol), 'color', currColor, 'LineWidth', 1);
                        plt.Color(4) = singleTrialAlpha;
                        meanAnnotData(startVol:endVol, iTrial, uniqueAnnotVals == currVolType) = currData(startVol:endVol); % --> [volume, trial, annotVal]   
                    end
                end
                
            else
                % Plot trial dF/F if not using annotation data for color coding
                if nGroups == 1
                    currColor = cm(iTrial, :);
                else
                    currColor = cm(iGroup, :);
                end
                plt = plot(ax, volTimes, smooth(currData, smoothWin), 'color', currColor, 'LineWidth', 1);
                plt.Color(4) = singleTrialAlpha;
                
            end%if
        end%iTrial
    end%if
    
    % Shade one SD above and below mean
    if shadeSDs
        upper = smooth(groupAvgDff, 3) + groupStdDev;
        lower = smooth(groupAvgDff, 3) - groupStdDev;
        jbfill(volTimes, upper', lower', shadeColor, shadeColor, 1, 0.2);
    end
    
    % Plot mean response line
    if ~isempty(annotArr)
        for iAnnot = 1:length(uniqueAnnotVals)
            currMeanAnnotData = mean(meanAnnotData(:,:,iAnnot), 2, 'omitnan'); % --> [volume]
            meanPlots(iAnnot) = plot(ax, volTimes, smooth(currMeanAnnotData, smoothWin), 'LineWidth', 2, 'Color', cm(uniqueAnnotVals(iAnnot) + 1, :));
        end
    end
    
    meanPlots(iGroup) = plot(ax, volTimes, smooth(groupAvgDff, smoothWin), 'LineWidth', 2, 'Color', meanLineColor * 1);
    
end%iGroup

% Plot alignment line if applicable
yL = ylim();
if volTimes(1) < 0
    plot(ax, [0 0], [-100 100], 'Color', 'k', 'LineWidth', 2); % Huge Y-range so it doesn't get cut off if I increase the ylims later
end
ylim(yL);

% Shade event timing if applicable
if ~isempty(eventShading)
    for iShade = 1:size(eventShading, 1)
        onsetTime = eventShading(iShade, 1);
        offsetTime = eventShading(iShade, 2);
        if size(eventShadeColor, 1) < size(eventShading, 1)
            shadeColor = eventShadeColor;
        else
            shadeColor = eventShadeColor(iShade, :);
        end
        plot_stim_shading([onsetTime, offsetTime], 'Axes', ax, 'Color', shadeColor);
        ylim(yL);
    end
end

% Add legend if applicable
if ~isempty(meanLineLegend)
    legend(meanPlots, meanLineLegend);
end

end