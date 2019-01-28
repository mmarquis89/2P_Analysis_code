function [plotHandle, plotAxes, plotFig] = plot_behavior_summary_2D(infoStruct, annotationArr, plotAxes, titleStr, trialGroups)
%=======================================================================================================
% PLOT A 2D SUMMARY FIGURE OF BEHAVIOR ANNOTATION DATA THROUGHOUT EXPERIMENT
% 
% Uses imagesc() to visualize annotated behavior data from an entire experiment. By default will plot
% trials in chronological order, but they can be optionally split into groups. Can create new figure 
% for the plot or plot in a specific axes object. 
%
% INPUTS:
%       infoStruct    = the main imaging data structure containing metadata for the experiment. 
%                         Specifically, must contain the fields "nFrames", "expDate", and "trialDuration".
%       annotationArr = array of behavioral annotation data to be plotted (row = trial, col = frame)
%       plotAxes      = <OPTIONAL> the handle to the axes you want the figure to be plotted in. Pass 
%                          [] to create a new figure and axes for the plot.
%       titleStr      = <OPTIONAL> string to be used as the title of the plot. Pass [] to use the default
%                         title of [experiment date, 'Behavior Summary']
%       trialGroups   = <OPTIONAL> nTrials x 1 numeric array specifying two or more trial groups (numbers
%                         starting at one) to split the data into. Higher numbered groups will be plotted
%                         below lower numbered groups. Pass [] to plot all trials in chronological order.
%       
% OUTPUTS:
%       plotHandle  = the handle to the plot that was created by the function
%       plotAxes    = the axes that the figure was plotted in
%       plotfig     = handle to the new figure, if one was created (otherwise returns [])
%========================================================================================================

% Remove trials from annotation array if necessary
goodTrials = infoStruct.goodTrials;
annotationArr(~goodTrials, :) = [];
if ~isempty(trialGroups)
    trialGroups(~goodTrials) = [];    
end



% Create or select figure and axes depending on whether an axes handle was provided
if isempty(plotAxes)
    plotFig = figure(1); clf; 
    plotAxes = axes();
    plotFig.Position = [200 100 1120 840];
    plotFig.Color = [1 1 1];
else
    axes(plotAxes)
    plotFig = [];
end

% Create custom colormap
cMap = [rgb('Indigo'); ...
        rgb('Orange'); ...
        rgb('Cyan'); ...
        rgb('Green'); ...
        rgb('Gold'); ...
        rgb('Yellow'); ...
        rgb('Green'); ...
        rgb('Red'); ...
        rgb('blue'); ...
        rgb('Gold'); ...
        rgb('black') ...
        ];
annotationArr(end+1, 1:11) = 1:11; % To keep the color mapping consistent

% Plot data
if isempty(trialGroups)
    % Plot all trials in order
    plotHandle = imagesc(annotationArr);
    yl = ylim();
    ylim(yl - [0 1]) % To cut off the extra row that's just there for colormapping
else
    % Separate by trial groups
    plotArr = [];
    spacerArr = ones(4, size(annotationArr, 2)) * (max(annotationArr(:)) - 1);
    for iGroup = 1:length(unique(trialGroups))
        plotArr = [plotArr; spacerArr; annotationArr(trialGroups == iGroup,:)];
    end
    plotHandle = imagesc(plotArr);
end
colormap(plotAxes, cMap)

% Format axes
plotAxes.FontSize = 12;
plotAxes.XLabel.String = 'Time (sec)';
plotAxes.YLabel.String = 'Trial number';
plotAxes.XLabel.FontSize = 16;
plotAxes.YLabel.FontSize = 14;
plotAxes.XTick = (0:1:infoStruct.trialDuration) * infoStruct.FRAME_RATE;
plotAxes.XTickLabel = 0:1:infoStruct.trialDuration;
if isempty(titleStr)
    plotAxes.Title.String = [regexprep(infoStruct.expDate, '_(?<num>..)', '\\_$<num>'), ...
                                'Behavior Summary']; % regex to add escape characters
else
    plotAxes.Title.String = titleStr;
end

end