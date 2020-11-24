%% Set directories and get expList
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData_60D05_7f';
figDir = fullfile(parentDir, 'Figs');

fileList = dir(fullfile(parentDir, '2020*'));
expList = unique(cellfun(@(x) x(1:10), {fileList.name}, 'uniformoutput', 0));

% ----- Load data for the selected experiments -----

% Load metadata 
[expMd, trialMd] = load_metadata(expList, parentDir);

% Load imaging data
roiData = load_roi_data(expList, parentDir);

% Load FicTrac data
ftData = load_ft_data(expList, parentDir);

% Load flailing event data (if any exists)
try
    eventData = load_event_data(expList, parentDir);
    flailingEvents = eventData.flailing;
catch; end

% Load panels metadata
panelsMetadata = load_panels_metadata(expList, parentDir);

% Consolidate wedge data and calculate bump population vector averages
try
    
% Calc dF/F for EB wedge ROIs
currData = roiData(~cellfun(@isempty, regexp(roiData.roiName, 'EB-', 'match')), :);
currData.trialDff = cellfun(@(x, y) (x - y) ./ y, currData.rawFl, ...
        mat2cell(currData.trialBaseline, ones(size(currData.trialBaseline))), 'uniformoutput', 0);
currData.expDff = cellfun(@(x, y) (x - y) ./ y, currData.rawFl, ...
        mat2cell(currData.expBaseline, ones(size(currData.expBaseline))), 'uniformoutput', 0);
    
wedgeData = [];
for iExp = 1:size(expMd, 1)
    currExpData = currData(strcmp(currData.expID, expMd.expID{iExp}), :);
    for iTrial = 1:max(currExpData.trialNum)
        
        % Consolidate Fl data into matrices
        td = sortrows(currExpData(currExpData.trialNum == iTrial, :), 'roiName');
        rawFlMat = cell2mat(td.rawFl');
        trialDffMat = cell2mat(td.trialDff');
        expDffMat = cell2mat(td.expDff');
        
        % Calculate bump population vector average
        angleMat = repmat((-7 * pi/8):pi/4:(7 * pi/8), size(expDffMat, 1), 1);
        [x, y] = pol2cart(angleMat, expDffMat); % Convert to cartesian coordinates to add vectors
        [theta, ~] = cart2pol(sum(x, 2), sum(y, 2));
        theta = -theta; % Inverting sign so it matches other data in plots 
        
        % Compile into table and append to wedgeData
        newRow = table(expMd.expID(iExp), iTrial, {rawFlMat}, {trialDffMat}, {expDffMat}, {theta}, ...
                {theta / pi * 4 + 4.5}, 'variableNames', {'expID', 'trialNum', 'rawFl', 'trialDff', ...
                'expDff', 'pvaRad', 'pvaWedge'});
        wedgeData = [wedgeData; newRow];
        
    end
end
catch ME; rethrow(ME); end

% Consolidate PB glomerulus data
try
    
% Calc dF/F for PB glom ROIs
currData = roiData(cellfun(@isempty, regexp(roiData.roiName, 'EB-', 'match')), :);
currData.trialDff = cellfun(@(x, y) (x - y) ./ y, currData.rawFl, ...
        mat2cell(currData.trialBaseline, ones(size(currData.trialBaseline))), 'uniformoutput', 0);
currData.expDff = cellfun(@(x, y) (x - y) ./ y, currData.rawFl, ...
        mat2cell(currData.expBaseline, ones(size(currData.expBaseline))), 'uniformoutput', 0);
    
glomData = [];
for iExp = 1:size(expMd, 1)
    currExpData = currData(strcmp(currData.expID, expMd.expID{iExp}), :);
    for iTrial = 1:max(currExpData.trialNum)
        
        % Consolidate Fl data into matrices
        td = sortrows(currExpData(currExpData.trialNum == iTrial, :), 'roiName');
        rawFlMat = cell2mat(td.rawFl');
        trialDffMat = cell2mat(td.trialDff');
        expDffMat = cell2mat(td.expDff');
        
        % Sort into paired order
        rawFlMat = rawFlMat(:, [1:9, 16:-1:10]);
        trialDffMat = trialDffMat(:, [1:9, 16:-1:10]);
        expDffMat = expDffMat(:, [1:9, 16:-1:10]);
        
        % Compile into table and append to glomData
        newRow = table(expMd.expID(iExp), iTrial, {rawFlMat}, {trialDffMat}, {expDffMat}, ...
                'variableNames', {'expID', 'trialNum', 'rawFl', 'trialDff', 'expDff'});
        glomData = [glomData; newRow];
        
    end
end
disp('Loading complete')
catch ME; rethrow(ME); end

%% CHECK CORRELATION BETWEEN GLOMERULI

expList = {'20201120-1', '20201120-2', '20201120-3'}%unique(expMd.expID);

try

glomPairNames = table((1:8)', {'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'}', ...
    {'R1', 'R8', 'R7', 'R6', 'R5', 'R4', 'R3', 'R2'}', 'variablenames', ...
    {'wedge', 'leftGlom', 'rightGlom'});

for iExp = 1:size(expMd, 1)
    
    
    currExpID = expMd.expID{iExp};
    
    if ismember(currExpID, expList)
        currExpRoiData = roiData(strcmp(roiData.expID, currExpID), :);
        
        % Concatenate data from all trials
        flMat = [];
        for iTrial = 1:max(currExpRoiData.trialNum)
            currTrialRoiData = currExpRoiData(currExpRoiData.trialNum == iTrial, :);
            
            currLeftFl = [];
            currRightFl = [];
            for iWedge = 1:8
                leftData = currTrialRoiData.rawFl(strcmp(currTrialRoiData.roiName, ...
                        glomPairNames.leftGlom{iWedge}));
                rightData = currTrialRoiData.rawFl(strcmp(currTrialRoiData.roiName, ...
                        glomPairNames.rightGlom{iWedge}));
                if ~isempty(leftData)
                    currLeftFl(:, iWedge) = leftData{:};
                else
                    currLeftFl(:, iWedge) = nan(size(currTrialRoiData.rawFl{1}));
                end
                if ~isempty(rightData)
                    currRightFl(:, iWedge) = rightData{:};
                else
                    currRightFl(:, iWedge) = nan(size(currTrialRoiData.rawFl{1}));
                end
            end
            currFl = [currLeftFl, currRightFl];
            flMat = [flMat; currFl];
        end
        
        % Calculate correlations between glomeruli
        R = corrcoef(flMat);
        LRCorrMat = R(1:8, 9:16);
        
        % Plot figure
        f = figure(iExp); clf;
        f.Color = [1 1 1];
        imagesc(LRCorrMat);
        ax = gca;
        colormap('bluewhitered')
        ax.YTick = 1:8;
        ax.XTick = 1:8;
        ax.YTickLabel = glomPairNames.leftGlom;
        ax.XTickLabel = glomPairNames.rightGlom;
        axis square
        title(currExpID)
    end
end
catch ME; rethrow(ME); end
    
   
%% 

expID = '20201117-1';

trialNum = 5;

sourceData = glomData;
sourceData = wedgeData;
plotPVA = 1;
plotMeanPVA = 1;
useFlow = 1;
flType = 'expDff';
flMax = [];
smWin = 5;
% trialDuration = 300;

currData = inner_join(sourceData, expMd, trialMd, panelsMetadata, ftData);

currExpData = currData(strcmp(currData.expID, expID), :);
td = currExpData(currExpData.trialNum == trialNum, :);

flMat = td.(flType){:};


% 50 = 344
% 65 = 447
% 90 = 619
% 130 = 894
% 195 = 1341

% flMat(1:447, :) = nan;
% flMat(448:894, :) = nan;
% flMat(895:1341, :) = nan;
% flMat(1342:end, :) = nan;

% flMat(1:619, :) = nan;
% flMat(620:812, :) = nan;
% flMat(813:1341, :) = nan;
% flMat(1342:end, :) = nan;


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
    flDataTemp = flMat;
    for iPos = 1:numel(unique(td.panelsPosX{:}))
        meanFlData(iPos, :) = ...
            mean(flDataTemp(panelsPosVols == (iPos - 1), :), 1, 'omitnan'); % --> [barPos, wedge]
    end
    meanFlShift = cat(1, meanFlData(92:96, :), meanFlData(1:91, :));
    cm = hsv(size(meanFlShift, 2)) .* 0.9;
end

f = figure(1); clf
f.Color = [1 1 1];
clear allAx
subaxis(5, 3, 1, 'mb', 0.05, 'mt', 0.03, 'ml', 0.05, 'mr', 0.03)

% Heatmap of mean visual tuning
allAx(1) = gca;
plotX = -180:3.75:(180 - 3.75);
if td.usingPanels
    smFl = smoothdata(meanFlShift, 1, 'gaussian', smWin);
    imagesc(plotX, [1 size(meanFlShift, 2)], smFl');
    hold on; plot([0 0], [0 size(meanFlShift, 2) + 1], 'color', 'k', 'linewidth', 1)
    allAx(1).XTick = -180:45:180;
    ylim([0.5 size(meanFlShift, 2) + 0.5])
    % xlabel('Bar position (deg)', 'fontsize', 10)
    ylabel('Visual tuning')
    
    % Overlay colored bars at edges to identify EB wedges across plots
    for iWedge = 1:size(cm, 1)
        plot(plotX(1:6), ones(1, 6) .* iWedge, 'color', cm(iWedge, :), 'linewidth', 4)
        plot(plotX(end - 5:end), ones(1, 6) .* iWedge, 'color', cm(iWedge, :), 'linewidth', 4)
    end
    
    % Overlay PVA for tuning
    if plotMeanPVA && size(flMat, 2) == 8
        angleMat = repmat((-7 * pi/8):pi/4:(7 * pi/8), size(smFl, 1), 1);
        [x, y] = pol2cart(angleMat, smFl); % Convert to cartesian coordinates to add vectors
        [theta, ~] = cart2pol(sum(x, 2), sum(y, 2));
        theta = -theta; % Inverting sign so it matches other data in plots
        theta = theta / pi * 4 + 4.5;
        plot(plotX, 9 - theta, 'color', 'g', 'linewidth', 2)
    end
    
end

% Plot of mean visual tuning
subaxis(5, 3, 2)
allAx(2) = gca;
hold on;
if td.usingPanels
    plotData = smoothdata(meanFlShift, 1, 'gaussian', smWin, 'omitnan');
    for iWedge = 1:size(meanFlShift, 2)
        plotX = -180:3.75:(180 - 3.75);
        plotX(isnan(meanFlShift(:, iWedge))) = nan;
        plot(plotX, plotData(:, iWedge), 'color', cm(iWedge, :), 'linewidth', 1);
    end
    allAx(2).XTick = -180:45:180;
    yL = [min(plotData(:)), max(plotData(:))] .* [0.9, 1.05];
    plot([0 0], yL, 'color', 'k', 'linewidth', 2)
    ylim(yL);
    xlim([-180, 180]);
end

% Mean visual tuning in polar coordinates
subaxis(5, 3, 3)
allAx(3) = gca;
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

% Mean dF/F data for each glomerulus
subaxis(5, 3, 4:6)
allAx(end + 1) = gca;
smMat = smoothdata(flMat, 1, 'gaussian', smWin);
smMat(1) = 0;
if ~isempty(flMax)
    smMat(end) = flMax;
    smMat(smMat > flMax) = flMax;
end
imagesc([0, td.trialDuration], [1, size(flMat, 2)], smMat')
hold on; 
% Overlay the dF/F population vector average
volTimes = (1:td.nVolumes) ./ td.volumeRate;
if size(flMat, 2) == 8
    if plotPVA
        plot(volTimes, 9 - smoothdata(td.pvaWedge{:}, 1, 'gaussian', smWin), 'color', 'g', ...
                'linewidth', 1.25)
    end
else
    plot([0, volTimes(end)], [8.5 8.5], 'color', 'g', 'linewidth', 3)
end
ylabel('dF/F and PVA');
colorbar; 
colormap(magma)

% Panels bar position
subaxis(5, 3, 7:9)
allAx(end + 1) = gca;
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
cb = colorbar; 
cb.Visible = 'off';

% FicTrac heading overlaid with dF/F pva
subaxis(5, 3, 10:12)
allAx(end + 1) = gca;
if useFlow
    currFlow = td.meanFlow{:};
    currFlow(end) = median(currFlow);
    currFlow(1:20) = median(currFlow);
    plotData = repeat_smooth(currFlow, 20, 'dim', 1, 'smwin', smWin);
%     plotData2 = repeat_smooth(currFlow, 20, 'dim', 2, 'smwin', 6);
    plotData = plotData - min(plotData);
%     plotData2 = plotData2 - min(plotData2);
    flowFrameDur = median(diff(td.frameTimes{:}));
    flowFrameTimes = (1:1:numel(currFlow)) * flowFrameDur;
    plot(flowFrameTimes, plotData, 'color', 'k', 'linewidth', 1);
    hold on;
%     plot(flowFrameTimes, plotData2, 'color', 'g');
    cb = colorbar;
    cb.Visible = 'off';
else
    HD = repeat_smooth(unwrap(td.intHD{:}), 20, 'smWin', smWin);
    plot(td.frameTimes{:}, 2*pi - mod(HD, 2*pi), 'color', 'k');
    ylim([0 2*pi])
    ylabel('Fly heading (rad)')
    yyaxis('right')
    uwVectAvgRad = smoothdata(unwrap(td.pvaRad{:} + pi), 1, 'gaussian',smWin);
    uwVectAvgRad = uwVectAvgRad - uwVectAvgRad(1) + 2*pi;
    plot(volTimes, mod(uwVectAvgRad, 2*pi), 'color', rgb('darkorange'));
    % ylabel('PVA')
    ylim([0 2*pi])
    xlim([0 volTimes(end)])
    cb = colorbar;
    cb.Visible = 'off';
    legend({'Heading', 'PVA'}, 'location', 'nw')
end

% FicTrac movement speed
subaxis(5, 3, 13:15)
allAx(end + 1) = gca;
plot(td.frameTimes{:}, repeat_smooth(td.moveSpeed{:}, 15, 'smWin', smWin), 'color', 'k');%colorbar
ylabel('Move speed (mm/sec)')
ylim([0 20])
yyaxis('right')
plot(td.frameTimes{:}, abs(repeat_smooth(td.yawSpeed{:}, 15, 'smWin', smWin)), ...
        'color', rgb('orange'));
allAx(end).YTick = [];
% ylabel('Yaw speed (rad/sec)');
cb = colorbar; 
cb.Visible = 'off';
xlabel('Time (s)');
legend({'Move speed', 'Yaw speed'}, 'location', 'nw')

% Link the X-axis limits across all plots
% linkaxes(allAx(3:end), 'x');
linkaxes(allAx([4 5 6 7]), 'x');
xlim([0 td.trialDuration])
















