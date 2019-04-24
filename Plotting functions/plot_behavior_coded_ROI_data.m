function ax = plot_behavior_coded_ROI_data(ax, flData, behavData, varargin)
%===================================================================================================
% PLOT MEAN dF/F DATA WITHIN ROI THROUGHOUT TRIAL COLORED BY BEHAVIOR
%
% INPUTS:
%
%       ax        = the handle to the axes where the plot should be created
%
%       flData    = the array of ROI-averaged fluorescence data with format: [volume, trial]
%
%       behavData = an array of behavior data (either annotation or FicTrac data) with size
%                   [m x nTrials]. If the size of the first dimension does not match flData, it
%                   will be resampled to match the volume rate, so the two arrays should represent
%                   the same period of time even if the sampling rates differ. The values will be
%                   used to color the corresponding flData values in the output plot. Datapoints
%                   with a value of nan will be omitted from the plot entirely.
%
% OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%
%       'OutlierSD'        = (default: 5) the number of standard deviations to use as an outlier
%                            exclusion threshold.
%
%       'SingleTrials'     = (default: 1) a boolean specifying whether to plot all individual trials
%                            in the background of the mean dF/F line. Should be set to one if 
%                            'EdgeColorMode' is set to 'flat'.
%
%       'StdDevShading'    = (default: 1) boolean specifying whether to shade 1 standard deviation
%                            above and below the mean flData line.
%
%
%       'SmoothWinSize'    = (default: 3) width in volumes of the moving average smoothing window
%                            that will be applied to BOTH the mean plot and the individual trial
%                            lines. Use a window size of "1" to skip smoothing entirely.
%
%       'VolOffset'        = (default: 0) a scalar number to add to the volumes for the purpose of
%                            labeling the X-axis. For example, if you are imaging at 10 volumes/sec
%                            and want your X-axis labeling to start at -2 seconds, use -20.
%
%       'SingleTrialAlpha' = (default: 1) a value from 0-1 specifying the transparency of the single
%                            trial lines.
%
%       'VolumeRate'       = (default: 6.87) the rate (in hz) at which imaging volumes were acquired
%
%       'Colormap'         = (default: 'parula') specifies a colormap to map the behavior data onto. 
%                            Can be either an Nx3 array of actual color values, or a valid colormap 
%                            name (e.g. 'parula').
%
%       'EdgeColorMode'    = (default: 'interp') specifies whether the behavior data is a set of 
%                            discrete categorical values (e.g. annotation data) or a continuous 
%                            variable (e.g. FicTrac data). Can be set to 'flat' for categorical
%                            data or 'interp' for continuous data. 
%===================================================================================================

% Parse optional arguments
p = inputParser;
addParameter(p, 'OutlierSD', 5);
addParameter(p, 'SingleTrials', 1);
addParameter(p, 'StdDevShading', 1);
addParameter(p, 'SmoothWinSize', 3);
addParameter(p, 'VolOffset', 0);
addParameter(p, 'SingleTrialAlpha', 1);
addParameter(p, 'VolumeRate', 6.87);
addParameter(p, 'Colormap', 'parula');
addParameter(p, 'EdgeColorMode', 'interp');
addParameter(p, 'YAxisLabel', 'dF/F');
parse(p, varargin{:});
outlierSD = p.Results.OutlierSD;
singleTrials = p.Results.SingleTrials;
shadeSDs = p.Results.StdDevShading;
smoothWin = p.Results.SmoothWinSize;
volumeRate = p.Results.VolumeRate;
volOffset = p.Results.VolOffset;
singleTrialAlpha = p.Results.SingleTrialAlpha;
cm = p.Results.Colormap;
yAxisLabel = p.Results.YAxisLabel;
edgeColorMode = p.Results.EdgeColorMode;

% --------------------------------------------------------------------------------------------------

% Setup variables
nVolumes = size(flData, 1);
nTrials = size(flData, 2);
volTimes = ([1:1:nVolumes] + volOffset) ./ volumeRate;

% Resample behavior data if necessary
if size(behavData, 1) ~= size(flData, 1)
    behavData = resample(behavData, size(flData, 1), size(behavData, 1)); % --> [volume, trial]
end

% Format axes
hold on
ax.YLabel.String =  yAxisLabel;
ax.XLabel.String = 'Time (s)';
xlim(ax, [min(volTimes), max(volTimes)]);

% Calculate average fl data
avgFl = mean(flData, 2, 'omitnan');

% Discard any trials that are too many SDs from mean
outliers = zeros(1, size(flData, 2));
for iTrial = 1:size(flData, 2)
        if max(abs(flData(:, iTrial) - avgFl) > (outlierSD * mean(std(flData, 0, 2))))
            outliers(iTrial) = 1;
        end
end
if sum(outliers) > 0
    disp(['Omitting ' num2str(sum(outliers)), ' outlier trials'])
end
flData(:, logical(outliers)) = [];      % --> [volume, trial]
behavData(:, logical(outliers)) = [];
stdDev = std(flData, 0, 2);

% Re-calculate average without outliers
avgFl = mean(flData, 2, 'omitnan');

% Plot individual trials
if singleTrials
    for iTrial = 1:size(flData, 2)
        currFl = movmean(flData(:, iTrial), smoothWin, 'omitnan');
        currBehavData = behavData(:, iTrial);
        currX = ((1:nVolumes) ./ volumeRate)';
        currY = currFl;
        currXY = [currX, currY];
        currXY = currXY(~any(isnan(currXY), 2), :);                 % to skip over any nan values
        currBehavData = currBehavData(~any(isnan(currXY), 2), :);   % to skip over any nan values
        surface('XData', [currXY(:,1) currXY(:,1)], ...
                'YData', [currXY(:,2) currXY(:,2)], ...
                'ZData', zeros(numel(currXY(:,1)), 2), ...
                'CData', [currBehavData, currBehavData] + 1, ...
                'FaceColor', 'none', ...
                'EdgeAlpha', singleTrialAlpha, ...
                'EdgeColor', edgeColorMode, ...
                'Marker', 'none');
    end
end

% Add colorbar if data is continuous
if strcmp(edgeColorMode, 'interp')
    colorbar(ax)
end

% Plot mean response line
if singleTrials
    plot(ax, volTimes, movmean(avgFl, smoothWin, 'omitnan'), 'LineWidth', 2, 'Color', 'k');
else
    avgBehavData = mean(behavData, 2, 'omitnan');
    avgBehavData = [avgBehavData(2); avgBehavData(2:end)]; % to drop artifically low first trial
    currX = ((1:nVolumes) ./ volumeRate)';
    currY = movmean(avgFl, smoothWin, 'omitnan');
    currXY = [currX, currY];
    currXY = currXY(~any(isnan(currXY), 2), :);             % to skip over any nan values
    avgBehavData = avgBehavData(~any(isnan(currXY), 2), :); % to skip over any nan values
    surface('XData', [currXY(:,1) currXY(:,1)], ...
            'YData', [currXY(:,2) currXY(:,2)], ...
            'ZData', zeros(numel(currXY(:,1)), 2), ...
            'CData', [avgBehavData, avgBehavData], ...
            'FaceColor', 'none', ...
            'EdgeColor', 'interp', ...
            'Marker', 'none', ...
            'LineWidth', 3);
end

% Set colormap
colormap(cm)

% If data is categorical, plot addtional mean lines for those volumes only
if strcmp(edgeColorMode, 'flat')
    behavVals = unique(behavData(~isnan(behavData)));
    for iVal = 1:numel(behavVals)
        currFlData = flData;
        currFlData(behavData ~= behavVals(iVal)) = nan;
        behavValMean = mean(currFlData, 2, 'omitnan');
        plot(ax, volTimes, movmean(behavValMean, smoothWin, 1, 'omitnan'), 'LineWidth', 2, 'Color', cm((behavVals(iVal)+1), :)*0.75)
    end
end

% Shade one SD above and below mean in grey
if shadeSDs
    upper = movmean(groupAvgFl, 3, 1, 'omitnan') + groupStdDev;
    lower = movmean(groupAvgFl, 3, 1, 'omitnan') - groupStdDev;
    jbfill(volTimes, upper', lower', shadeColor, shadeColor, 1, 0.2);
end

% Plot alignment line if applicable
yL = ylim();
if volTimes(1) < 0
    plot(ax, [0 0], [-100 100], 'Color', 'k', 'LineWidth', 2); % Huge Y-range so it doesn't get cut off if I increase the ylims later
end
ylim(yL);


end