function [f, allAx] = plot_single_trial_visual_tuning_summary(trialData, plotParams)
%===================================================================================================
% Plots a large summary figure of the visual tuning of EB wedges or PB glomeruli throughout a trial 
% of open loop bar swinging. Includes an averaged tuning heatmap, overlaid tuning curves both 
% cartesian and polar plots, and the raw imagesc plot of imaging data with aligned plots of bar 
% location and fly movement variables (for figure sizing consistency, the bottom subplot will be 
% left blank for experiments in which the fly was not mounted on the ball). 
% 
% INPUTS:
% 
%       trialData = a table consisting of the joined tables of exp/trial metadata, panels metadata, 
%                   FicTrac data, and either wedge or glomerulus imaging data as created by 
%                   the load_PB_data() function. The specific fields of the table that will be used 
%                   by this function are as follows:
%                           rawFl, trialDff, or expDff (depending on p.flType)
%                           pvaRad       (if p.plotPVA)
%                           pvaWedge     (if p.plotPVA)
%                           nVolumes 
%                           volumeRate
%                           trialDuration
%                           panelsDisplayRate
%                           nPanelsFrames
%                           usingPanels
%                           panelsPosX   (if usingPanels)
%                           intHD        (if ~p.useFlow)
%                           moveSpeed
%                           yawSpeed
%                           frameTimes
%                           vidFrameTimes(if p.useFlow)
%                           meanFlow     (if p.useFlow)
%       
%       plotParams = a struct containing the following fields of plotting parameters:
%                           flType      ('rawFl', 'trialDff', or 'expDff')
%                           plotPVA     (boolean, overlay PVA on main imagesc plot?)
%                           plotMeanPVA (boolean)
%                           useFlow     (boolean)
%                           flMax       (value to cap fl at for imagesc plots, or [] to skip)
%                           smWin       (width of gaussian smoothing window)
%                           figNum      (number of figure to create, or [] to default to 1)
% 
% OUTPUTS:
% 
%       f       = handle of the figure that was created               
% 
%       allAx   = vector of handles the subplot axes handles
% 
%===================================================================================================

td = trialData;
p = plotParams;
flMat = td.(p.flType){:};

if isempty(p.figNum)
   p.figNum = 1; 
end

% Create figure
f = figure(p.figNum); clf
f.Color = [1 1 1];
subaxis(5, 3, 1, 'mb', 0.06, 'mt', 0.03, 'ml', 0.05, 'mr', 0.03);
    
% Get mean panels pos data
panelsFrameTimes = double(1:td.nPanelsFrames) ./ td.panelsDisplayRate;
volTimes = (1:td.nVolumes) ./ td.volumeRate;
panelsPosVols = [];
if td.usingPanels
    for iVol = 1:size(flMat, 1)
        [~, currVol] = min(abs(panelsFrameTimes - volTimes(iVol)));
        panelsPosVols(iVol) = td.panelsPosX{:}(currVol);
    end
    meanFlData = [];
    for iPos = 1:numel(unique(td.panelsPosX{:}))
        meanFlData(iPos, :) = ...
            mean(flMat(panelsPosVols == (iPos - 1), :), 1, 'omitnan'); % --> [barPos, wedge]
    end
    % Shift Fl data so that the center index is directly in front of the fly
    meanFlShift = cat(1, meanFlData(92:96, :), meanFlData(1:91, :));     
    cm = hsv(size(meanFlShift, 2)) .* 0.9;
    
    %  ------------- Heatmap of mean visual tuning -------------------------------------------------
    allAx(1) = gca();
    plotX = -180:3.75:(180 - 3.75);
    smFl = smoothdata(meanFlShift, 1, 'gaussian', p.smWin);
    imagesc(plotX, [1 size(meanFlShift, 2)], smFl');
    hold on; plot([0 0], [0 size(meanFlShift, 2) + 1], 'color', 'k', 'linewidth', 1)
    allAx(1).XTick = -180:45:180;
    ylim([0.5 size(meanFlShift, 2) + 0.5])
    % xlabel('Bar position (deg)', 'fontsize', 10)
    ylabel('Visual tuning')
    allAx(1).FontSize = 12;

    % Overlay colored bars at edges to identify EB wedges across plots
    for iWedge = 1:size(cm, 1)
        plot(plotX(1:6), ones(1, 6) .* iWedge, 'color', cm(iWedge, :), 'linewidth', 4)
        plot(plotX(end - 5:end), ones(1, 6) .* iWedge, 'color', cm(iWedge, :), 'linewidth', 4)
    end

    % Overlay PVA for tuning
    if p.plotMeanPVA && size(flMat, 2) == 8
        angleMat = repmat((-7 * pi/8):pi/4:(7 * pi/8), size(smFl, 1), 1);
        [x, y] = pol2cart(angleMat, smFl); % Convert to cartesian coordinates to add vectors
        [theta, ~] = cart2pol(sum(x, 2), sum(y, 2));
        theta = -theta; % Inverting sign so it matches other data in plots
        theta = theta / pi * 4 + 4.5;
        plot(plotX, 9 - theta, 'color', 'g', 'linewidth', 2)
    end
    
    %  ------------- Plot of mean visual tuning ----------------------------------------------------
    allAx(2) = subaxis(5, 3, 2);
    hold on;
    plotData = smoothdata(meanFlShift, 1, 'gaussian', p.smWin, 'omitnan');
    for iWedge = 1:size(meanFlShift, 2)
        plotX = -180:3.75:(180 - 3.75);
        plotX(isnan(meanFlShift(:, iWedge))) = nan;
        plot(plotX, plotData(:, iWedge), 'color', cm(iWedge, :), 'linewidth', 1);
    end
    allAx(2).XTick = -180:45:180;
    allAx(2).YTick = [];
    yL = [min(plotData(:)), max(plotData(:))] .* [0.9, 1.05];
    plot([0 0], yL, 'color', 'k', 'linewidth', 2)
    ylim(yL);
    xlim([-180, 180]);
    allAx(2).FontSize = 12;
    
end%if

% ------------- Mean visual tuning in polar coordinates --------------------------------------------
allAx(3) = subaxis(5, 3, 3);
plotX = deg2rad(-180:3.75:(180 - 3.75));
pax = polaraxes(); hold on
pax.Position = allAx(3).Position;
pax.Position = pax.Position - [0.04 0.025 0 0];
pax.Position = pax.Position .* [1 1 1.2 1.2];
allAx(3).Visible = 'off';
if td.usingPanels
    for iWedge = 1:size(meanFlShift, 2)
        plotData = smoothdata(meanFlShift(:, iWedge), 1, 'gaussian', 6) ;
        plotData = plotData - min(meanFlShift(:));
        polarplot(plotX, plotData, 'color', cm(iWedge, :), 'linewidth', 1);
    end
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.ThetaTick = [0:45:179, 225:45:359];
    pax.ThetaTickLabel = {'0' '+45', '+90' '+135' '-135', '-90', '-45'};
    pax.RTick = [];
    pax.FontSize = 12;
end

% ------------- Mean dF/F data for each glomerulus -------------------------------------------------
allAx(4) = subaxis(5, 3, 4:6);
smMat = smoothdata(flMat, 1, 'gaussian', p.smWin);
smMat(1) = 0;
if ~isempty(p.flMax)
    smMat(end) = p.flMax;
    smMat(smMat > p.flMax) = p.flMax;
end
imagesc([0, td.trialDuration], [1, size(flMat, 2)], smMat')
hold on; 

% Overlay colored bars at edges to identify EB wedges across plots
barLen = round(td.nVolumes ./ 70);
for iWedge = 1:size(cm, 1)
    plot(volTimes(1:barLen), ones(1, barLen) .* iWedge, 'color', cm(iWedge, :), 'linewidth', 4)
    plot(volTimes(end - (barLen-1):end), ones(1, barLen) .* iWedge, 'color', cm(iWedge, :), ...
            'linewidth', 4)
end

% Overlay the dF/F population vector average
volTimes = (1:td.nVolumes) ./ td.volumeRate;
if size(flMat, 2) == 8
    if p.plotPVA
        % Don't plot the beginning or end of the PVA so as to not overlap with the colored bars at 
        % the edges of the plot
        plot(volTimes((barLen+1):end-barLen), 9 - smoothdata(td.pvaWedge{:}((barLen+1):end-barLen), ...
                1, 'gaussian', p.smWin), 'color', 'g', 'linewidth', 1.25)
    end
else
    plot([0, volTimes(end)], [8.5 8.5], 'color', 'g', 'linewidth', 3)
end
    
ylabel('dF/F and PVA');
allAx(4).FontSize = 12;
colorbar; 
colormap(magma)

% ------------- Panels bar position ----------------------------------------------------------------
allAx(5) = subaxis(5, 3, 7:9);
if td.usingPanels
    plot(panelsFrameTimes, td.panelsPosX{:}, 'linewidth', 1.5, 'color', 'b');%colorbar
    barCenteredFrames = find(td.panelsPosX{:} == 44);
    yL = ylim;
    hold on;
    plotX = panelsFrameTimes(barCenteredFrames);
    plot([plotX; plotX], repmat(yL', 1, numel(barCenteredFrames)), 'color', ...
            'r');
end
ylabel('Bar position')
allAx(end).YTickLabel = [];
allAx(5).FontSize = 12;
cb = colorbar; 
cb.Visible = 'off';

% ------------- FicTrac heading overlaid with dF/F pva ---------------------------------------------
allAx(6) = subaxis(5, 3, 10:12);
if p.useFlow
    currFlow = td.meanFlow{:};
    if ~isempty(currFlow)
        currFlow(end) = median(currFlow);
        currFlow(1:20) = median(currFlow);
        plotData = repeat_smooth(currFlow, 20, 'dim', 1, 'smwin', p.smWin);
        %     plotData2 = repeat_smooth(currFlow, 20, 'dim', 2, 'smwin', 6);
        plotData = plotData - min(plotData);
        %     plotData2 = plotData2 - min(plotData2);
        plot(td.vidFrameTimes{:}, plotData, 'color', 'k', 'linewidth', 1);
        hold on;
        %     plot(td.vidFrameTimes{:}, plotData2, 'color', 'g');
        cb = colorbar;
        cb.Visible = 'off';
    end
else
    HD = repeat_smooth(unwrap(td.intHD{:}), 20, 'smWin', p.smWin);
    plot(td.frameTimes{:}, 2*pi - mod(HD, 2*pi), 'color', 'k');
    ylim([0 2*pi])
    ylabel('Fly heading (rad)')
    yyaxis('right')
    if p.plotPVA
        uwVectAvgRad = smoothdata(unwrap(td.pvaRad{:} + pi), 1, 'gaussian',p.smWin);
        uwVectAvgRad = uwVectAvgRad - uwVectAvgRad(1) + 2*pi;
        plot(volTimes, mod(uwVectAvgRad, 2*pi), 'color', rgb('darkorange'));
    end
    % ylabel('PVA')
    ylim([0 2*pi])
    xlim([0 volTimes(end)])
    cb = colorbar;
    cb.Visible = 'off';
    legend({'Heading', 'PVA'}, 'location', 'nw')
end
allAx(6).FontSize = 12;

% ------------- FicTrac movement speed -------------------------------------------------------------
if ~p.useFlow
    allAx(7) = subaxis(5, 3, 13:15);
    plot(td.frameTimes{:}, repeat_smooth(td.moveSpeed{:}, 15, 'smWin', p.smWin), 'color', 'k');
    ylabel('Move speed (mm/sec)')
    ylim([0 20])
    yyaxis('right')
    plot(td.frameTimes{:}, abs(repeat_smooth(td.yawSpeed{:}, 15, 'smWin', p.smWin)), ...
            'color', rgb('orange'));
    allAx(end).YTick = [];
    % ylabel('Yaw speed (rad/sec)');
    cb = colorbar;
    cb.Visible = 'off';
    xlabel('Time (s)');
    legend({'Move speed', 'Yaw speed'}, 'location', 'nw')
    
    % Link the X-axis limits across time series plots
    linkaxes(allAx([4:7]), 'x');
    xlim([0 td.trialDuration])
    allAx(7).FontSize = 12;
else
    % Link the X-axis limits across time series plots
    linkaxes(allAx([4:6]), 'x');
    xlim([0 td.trialDuration])
end
xlabel('Time (sec)');
end%function