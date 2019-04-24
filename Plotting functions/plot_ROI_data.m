function ax = plot_ROI_data(ax, flData, varargin)
%===================================================================================================
% PLOT MEAN dF/F DATA WITHIN ROI THROUGHOUT TRIAL
%
% Creates a plot of trial-averaged ROI fluorescence data for an entire trial in a specified set of 
% axes. There are various optional additions such as single trial plotting and standard deviation 
% shading. Data is lightly smoothed by default but this can be disabled.
%
% INPUTS:
%
%       ax       = the handle to the axes where the plot should be created
%
%       flData   = the array of ROI-averaged fluorescence data with format: [volume, trial]
%
%
% OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%
%       'TrialGroups'      = (default: []) 1 x nTrials numeric array specifying two or more trial 
%                            groups (numbered starting at one) to split the data into. Trials 
%                            assigned to group zero will be omitted. Each group will be assigned
%                            a color (individual trials will no longer be plotted chronologically),
%                            and instead of plotting the average of all trials each group will
%                            have its own trial-averaged line.
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
%       'SmoothWinSize'    = (default: 3) width in volumes of the moving average smoothing window
%                            to be applied to BOTH the mean dF/F plot and the individual trial
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
%       'VolumeRate'       = (default: 6.87) rate (in hz) at which imaging volumes were acquired
%
%       'YAxisLabel'       = (default: 'dF/F') text to label the plot's Y-axis
%
%       'Colormap'         = alternative colormap to use for plotting trial groups
%
%===================================================================================================

% Parse optional arguments
p = inputParser;
addParameter(p, 'TrialGroups', []);
addParameter(p, 'OutlierSD', 5);
addParameter(p, 'SingleTrials', 1);
addParameter(p, 'StdDevShading', 1);
addParameter(p, 'SmoothWinSize', 3);
addParameter(p, 'VolOffset', 0);
addParameter(p, 'SingleTrialAlpha', 1);
addParameter(p, 'Legend', []);
addParameter(p, 'VolumeRate', 6.87);
addParameter(p, 'Colormap', []);
addParameter(p, 'YAxisLabel', 'dF/F');
parse(p, varargin{:});
trialGroups = p.Results.TrialGroups;
outlierSD = p.Results.OutlierSD;
singleTrials = p.Results.SingleTrials;
shadeSDs = p.Results.StdDevShading;
smoothWin = p.Results.SmoothWinSize;
volumeRate = p.Results.VolumeRate;
volOffset = p.Results.VolOffset;
singleTrialAlpha = p.Results.SingleTrialAlpha;
meanLineLegend = p.Results.Legend;
cm = p.Results.Colormap;
yAxisLabel = p.Results.YAxisLabel;

% --------------------------------------------------------------------------------------------------

% Setup variables
nVolumes = size(flData, 1);
volTimes = ([1:1:nVolumes] + volOffset) ./ volumeRate;
nTrials = size(flData, 2);
if ~isempty(trialGroups)
    nGroups = numel(unique(trialGroups(trialGroups ~= 0)));
else
    nGroups = 1;
    trialGroups = ones(1, nTrials);
end

% Format axes
hold on
ax.YLabel.String = yAxisLabel;
ax.XLabel.String = 'Time (s)';
xlim(ax, [min(volTimes), max(volTimes)]);

% Create colormap
if isempty(cm)
    if nGroups == 1
        cm = jet(nTrials);
    else
        cm = [rgb('blue'); ... 
                rgb('red'); ... 
                rgb('green'); ... 
                rgb('magenta'); ... 
                rgb('cyan'); ... 
                rgb('gold'); ...
                rgb('DarkRed'); ...
                rgb('Yellow'); ...
                rgb('lime'); ...
                rgb('black'); ...
                rgb('maroon'); ...
                rgb('grey')];
    end
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
    
    % Separate data from current trial group and calculate average fl data 
    groupAvgFl = mean(flData(:, trialGroups == iGroup), 2, 'omitnan');
    groupFlData = flData(:, trialGroups == iGroup); % --> [volume, trial]
    
    % Discard any trials that are too many SDs from mean
    outliers = zeros(1, size(flData, 2));
    for iTrial = 1:size(flData, 2)
        if trialGroups(iTrial) == iGroup
            if max(abs(flData(:, iTrial) - groupAvgFl) > (outlierSD * mean(std(flData, 0, 2, 'omitnan'))))
                outliers(iTrial) = 1;
            end
        end
    end
    groupOutliers = outliers(trialGroups == iGroup);
    if sum(groupOutliers) > 0
        disp(['Omitting ' num2str(sum(groupOutliers)), ' outlier trials'])
    end
    groupFlData(:, logical(groupOutliers)) = []; % --> [volume, trial]
    groupStdDev = std(groupFlData, 0, 2);
    
    % Re-calculate average without outliers
    groupAvgFl = mean(flData(:,  logical((trialGroups == iGroup) .* ~logical(outliers))), 2, 'omitnan');
    
    % Plot individual trials in background
    if singleTrials
        for iTrial = 1:size(groupFlData, 2)
            currData = movmean(groupFlData(:, iTrial), smoothWin, 1, 'omitnan');
            
            % Plot trial dF/F if not using annotation data for color coding
            if nGroups == 1
                currColor = cm(iTrial, :);
            else
                currColor = cm(iGroup, :);
            end
            plt = plot(ax, volTimes, movmean(currData, smoothWin, 1, 'omitnan'), 'color', currColor, 'LineWidth', 1);
            plt.Color(4) = singleTrialAlpha;
            
        end%iTrial
    end%if

    if shadeSDs
        upper = movmean(groupAvgFl, 3, 1, 'omitnan') + groupStdDev;
        lower = movmean(groupAvgFl, 3, 1, 'omitnan') - groupStdDev;
        jbfill(volTimes, upper', lower', shadeColor, shadeColor, 1, 0.2);
    end
    
    % Plot mean response line    
    meanPlots(iGroup) = plot(ax, volTimes, movmean(groupAvgFl, smoothWin, 'omitnan'), 'LineWidth', 2, 'Color', meanLineColor * 1);

end%iGroup

% Plot alignment line if applicable
yL = ylim();
if volTimes(1) < 0
    plot(ax, [0 0], [-100 100], 'Color', 'k', 'LineWidth', 2); % Huge Y-range so it doesn't get cut off if I increase the ylims later
end
ylim(yL);

% Add legend if applicable
if ~isempty(meanLineLegend)
    legend(meanPlots, meanLineLegend, 'autoupdate', 'off');
end

end