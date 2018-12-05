

%% FICTRAC VARIABLE PLOTTING SCRIPT

currTrial = 18;
% currData = selectData(:, [1 6 2 3 5 4],currTrial); % re-order columns to be [Frame Count, Seq num, xPos, yPos, Speed, HD]
currAllData = allData(:, :, currTrial); % --> [frame, variable]

% Split into individual components
frameCounter = squeeze(allData(:, 1, :));
dRotCam = squeeze(allData(:, [2 3 4], :));          % [X, Y, Z]
dRotError = squeeze(allData(:, 5, :));
dRotLab = squeeze(allData(:, [6 7 8], :));          % [X, Y, Z]
absOrientCam = squeeze(allData(:, [9 10 11], :));   % [X, Y, Z]
absOrientLab = squeeze(allData(:, [12 13 14], :));  % [X, Y, Z]
mmXY = squeeze(allData(:, [15 16], :));            % [X, Y]
intHD = squeeze(allData(:, 17, :));
moveDirLab = squeeze(allData(:, 18, :));
moveSpeed = squeeze(allData(:, 19, :));
intForwardMove = squeeze(allData(:, 20, :));
intSideMove = squeeze(allData(:, 21, :));
timestamp = squeeze(allData(:, 22, :));
seqNum = squeeze(allData(:, 23, :));

% figure(1);clf;hold on; ax = gca();
% plot(dRotCam(:,:,currTrial))
% legend 1 2 3
% title dRotCam
% xlabel('Time (sec)')
% ax.XTick = xTickFR(1:2:end);
% ax.XTickLabel = xTickLabels(1:2:end);
% xlim([0 size(dRotCam, 1)])
% % tightfig
% 
% figure(3);clf;hold on; ax = gca()
% plot(dRotLab(:,:,currTrial))
% % plot(filtfilt(kb, ka, dRotLab))
% legend 1 2 3
% title dRotLab
% ax.XTick = xTickFR;
% ax.XTickLabel = xTickLabels;
% % tightfig
%
% figure(4);clf;hold on; ax = gca();
% plot(absOrientCam(:,:,currTrial))
% legend 1  2 3 
% title absOrientCam
% ax.XTick = xTickFR(1:2:end);
% ax.XTickLabel = xTickLabels(1:2:end);
% xlabel('Time (sec)')
% xlim([0 size(absOrientCam, 1)])
% % tightfig
%
% figure(5);clf;hold on; ax = gca()
% plot(absOrientLab(:,:,currTrial))
% legend 1 2 3
% title absOrientLab
% ax.XTick = xTickFR;
% ax.XTickLabel = xTickLabels;
% % tightfig

% figure(6);clf;hold on; ax = gca();
% plot(dRotError(:,currTrial))
% legend X Y x y
% title dRotError
% ax.XTick = xTickFR;
% ax.XTickLabel = xTickLabels;
% % tightfig
% 
% figure(7);clf;hold on; ax = gca();
% uwHD = unwrap(intHD(:,currTrial));
% smHD = smooth(uwHD, 1);
% plot(smHD) %(mod(smHD, (2*pi)))
% legend HD hd
% title intHD
% ax.XTick = xTickFR;
% ax.XTickLabel = xTickLabels;
% % tightfig
%
figure(9);clf;hold on
plot(filtfilt(kb, ka, smooth(moveSpeed(:,currTrial), 7)))
title moveSpeed
% tightfig

% figure(10);clf;hold on; ax = gca();
% plot([intForwardMove(:,currTrial), intSideMove(:,currTrial), intHD(:,currTrial)])
% legend Forward Side HD
% title intForward+SideMove
% ax.XTick = xTickFR;
% ax.XTickLabel = xTickLabels;
% tightfig

% Movement map
figure(11); clf; hold on;
nSegs = 20;
vectorLen = floor(size(mmXY(:,:,currTrial), 1) / nSegs);
cm = jet(nSegs);
for iSeg = 1:nSegs
    if iSeg == 1
        pad = 1;
    else
        pad = 0;
    end
    currSegInd = (pad + ((iSeg-1)*vectorLen)):(iSeg * vectorLen);
    plot(smooth(mmXY(currSegInd, 1, currTrial),3), smooth(mmXY(currSegInd, 2, currTrial),3), 'Color', cm(iSeg, :))
end
lims = 1.1 * max(max(abs(mmXY(:,:,currTrial))));
xlim([-lims lims])
ylim([-lims lims])
axis image

% % Anvil annotations
% figure(12); clf; hold on
% plot(behaviorAnnotArr(currTrial, :), '-*', 'Color', 'k')
% cm = [rgb('Navy'); rgb('Cyan'); rgb('maroon'); rgb('Gold')];
% for iPlot = 1:4
%     currData = behaviorAnnotArr(currTrial, :);
%     currData(currData ~= iPlot - 1) = nan;
%     plot(currData, '*', 'color', cm(iPlot, :))
% end
% ylim([-1 4])
% ax = gca;
% ax.YTick = [0 1 2 3];
% ax.YTickLabel = {'Quiescence', 'Locomotion', 'Grooming', 'IsolatedMovement'};

%% CREATE FICTRAC + BEHAVIOR VIDS
vidDir = 'B:\Dropbox (HMS)\2P Data\Behavior Vids\2018_06_29\_Movies'

vidFiles = dir(fullfile(vidDir, ['sid*tid*.mp4']));

for iTrial = 1:n
    
    vidName = vidFiles(iTrial).name;
    disp(vidName)
    vidFile = fullfile(vidDir, vidName);
    currData = selectData(:, [1 6 2 3 5 4], iTrial); % --> [Frame Count, Seq num, xPos, yPos, Speed, HD]
    
    % Convert xy data from radians to mm
    r = 4.5;
    mmData = [currData(:,1:2), currData(:, 3:5) * r, currData(:, 6)];
    
    % Smooth velocity data
    smoothVelData = [];
    for iAxis = 1:size(mmData, 2)
        smoothWin = 3;
        smoothVelData(:, iAxis) = smooth(mmData(:, iAxis), smoothWin);  % --> [sample, axis, trial]
    end
    % LP filter velocity data
    rate = 2 * (11/25);
    [kb, ka] = butter(2,rate);
    velData = filtfilt(kb, ka, smoothVelData);
    
    % Heading direction data
    HD = currData(:,end);
    
    % Create VidReader and VidWriter
    myVid = VideoReader(vidFile);
    myVidWriter = VideoWriter(fullfile(parentDir, ['Combined_', vidName]), 'MPEG-4');
    myVidWriter.FrameRate = FRAME_RATE;
    open(myVidWriter)
    
    h = figure(10);
    for iFrame = 1:nFrames
        
        currFrame = uint8(readFrame(myVid));
        %
        clf
        ySize = size(currFrame, 1);
        xSize = size(currFrame, 2);
        
        % Create figure
        screenSize = get( groot, 'Screensize' );
        h.Color = [1 1 1];
        h.OuterPosition = [50 100 (screenSize(3) - 100) screenSize(4) - 150];
        xyRatio = h.Position(3) / h.Position(4);
        
        % Create axes
        clf
        M = 0.006;
        P = 0.00;
        axVid = subaxis(3,6,[1 2 7 8], 'S', 0, 'M', M, 'PB', 0.05);
        axMove = subaxis(3,6,[3 4 9 10], 'S', 0, 'M', M, 'PB', 0.06, 'PR', 0.01);
        axHD = subaxis(3,6,[5 6 11 12], 'S', 0, 'M', M, 'PB', 0.06, 'PL', 0.01);
        axYawSpeed = subaxis(3,6,[13:15], 'S', 0, 'M', M, 'PB', 0.05, 'PR', 0.01, 'PL', 0.008);
        axVel = subaxis(3,6,[16:18], 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.02);
        
        % Movie frame plot
        imshow(currFrame, [], 'Parent', axVid);
        axis off
        axVid.XTickLabel = [];
        
        % Movement map plot
        axes(axMove);
        hold on
        nSegs = trialDuration;
        vectorLen = floor(size(mmData, 1) / nSegs);
        cm = jet(nSegs);
        for iSeg = 1:nSegs
            if iSeg == 1
                pad = 1;
            else
                pad = 0;
            end
            currSegInd = (pad + ((iSeg-1)*vectorLen)):(iSeg * vectorLen);
            plot(axMove, mmData(currSegInd, 3), mmData(currSegInd, 4), 'Color', cm(iSeg, :))
        end
        axis equal; %axis square
        lims = 1.1 * max(abs([ylim(axMove), xlim(axMove)]));
        if lims < 1
            lims = 1;
        end
        xlim([-lims lims])
        ylim([-lims lims])
        
        legend({'2D movement (mm)'}, 'FontSize', 11)
        x = mmData(iFrame, 3);
        y = mmData(iFrame, 4);
        [arrowVec(1), arrowVec(2)] = pol2cart(HD(iFrame)+(pi * 1.5), 0.5);
        %         [testVec(1), testVec(2)] = pol2cart(test(iFrame) + (pi * 1.5), 0.5);
        arrow([x - arrowVec(1)/2, y-arrowVec(2)/2], [x + arrowVec(1)/2, y + arrowVec(2)/2], 'FaceColor', 'green');
        %         arrow([x - testVec(1)/2, y-testVec(2)/2], [x + testVec(1)/2, y + testVec(2)/2], 'FaceColor', 'red');
        
        % Head direction plot
        axes(axHD)
        plt = polarplot([HD(iFrame), HD(iFrame)], [0 1], '-', 'LineWidth', 5, 'Color', 'k');
        axHD = plt.Parent;
        axHD.FontSize = 12;
        axHD.ThetaZeroLocation = 'bottom';
        axHD.RTickLabel = [];
        axHD.RTick = [];
        axes(axHD);
        hold on
        polarplot(linspace(0, pi * 2, 100), ones(100,1), 'linewidth', 0.5, 'Color', 'k');
        
        % Displacement plot
        axes(axYawSpeed);
        hold on
        plot(mmData(:, [3 4]));
        axYawSpeed.XTick = xTickFR;
        axYawSpeed.XTickLabel = xTickLabels;
        legend({'X disp (mm)', 'Y disp (mm)'}, 'FontSize', 11, 'Location', 'NorthWest')
        xlabel('Time (sec)')
        plot(iFrame, mmData(iFrame, 3), 'p', 'Color', 'k', 'MarkerSize', 12, 'MarkerFaceColor', 'k')
        plot(iFrame, mmData(iFrame, 4), 'p', 'Color', 'k', 'MarkerSize', 12, 'MarkerFaceColor', 'k')
        
        % Velocity plot
        axes(axVel)
        plot(velData(:,5));
        axVel.XTick = xTickFR;
        axVel.XTickLabel = xTickLabels;
        legend({'Speed (mm/sec)'}, 'FontSize', 11)
        xlabel('Time (sec)')
        hold on
        plot(iFrame, velData(iFrame, 5), 'p', 'Color', 'k', 'MarkerSize', 12, 'MarkerFaceColor', 'k')
        
        % Write frame to video(s)
        writeFrame = getframe(h);
        writeVideo(myVidWriter, writeFrame);
    end%iFrame
    
    close(myVidWriter)
end


%% Calibration video data

allData = csvread('calibration.dat');

frameCounter = squeeze(allData(:, 1));
dRotCam = squeeze(allData(:, [2 3 4]));          % [X, Y, Z]
dRotError = squeeze(allData(:, 5));
dRotLab = squeeze(allData(:, [6 7 8]));          % [X, Y, Z]
absOrientCam = squeeze(allData(:, [9 10 11]));   % [X, Y, Z]
absOrientLab = squeeze(allData(:, [12 13 14]));  % [X, Y, Z]
mmXY = squeeze(allData(:, [15 16]));            % [X, Y]
intHD = squeeze(allData(:, 17));
moveDirLab = squeeze(allData(:, 18));
moveSpeed = squeeze(allData(:, 19));
intFWSideMove = squeeze(allData(:, [20 21]));
intForwardMove = squeeze(allData(:, 20));
intSideMove = squeeze(allData(:, 21));
timestamp = squeeze(allData(:, 22));
seqNum = squeeze(allData(:, 23)); 

vidDir = '..'
vidFile = dir(fullfile(vidDir, 'calibration.avi'));
vidName = vidFile.name;
vidFile = fullfile(vidDir, vidName);

% Invert heading direction data
HD = (2*pi) - intHD;

% Create VidReader and VidWriter
myVid = VideoReader(vidFile);
myVidWriter = VideoWriter(fullfile(vidDir, ['Combined_', vidName, '_2']), 'MPEG-4');
myVidWriter.FrameRate = FRAME_RATE;
open(myVidWriter)

% Create figure
h = figure(100);
screenSize = get( groot, 'Screensize' );
h.Color = [1 1 1];
h.OuterPosition = [50 100 (screenSize(3) - 100) screenSize(4) - 150];

for iFrame = 1:1863
    
    currFrame = uint8(readFrame(myVid));
    %
    clf
    ySize = size(currFrame, 1);
    xSize = size(currFrame, 2);
    

    
    % Create axes
    clf
    M = 0.001;
    P = 0.01;
    axG = [6 12];
    markerSize = 10;
    
    axVid =             subaxis(axG(1), axG(2), 1,1, 3,3, 'S', 0, 'M', M);
    ax2DPos =           subaxis(axG(1), axG(2), 1,4, 3,3, 'S', 0, 'M', M, 'PB', 0.03, 'PL', 0.02);
    
    axMoveDirLab =      subaxis(axG(1), axG(2), 4,2, 1,2, 'S', 0, 'M', M, 'PT', 0.01, 'PR', 0.01);
    axHD =              subaxis(axG(1), axG(2), 4,4, 1,2, 'S', 0, 'M', M, 'PT', 0.01, 'PR', 0.01);   

    axdRotCam =         subaxis(axG(1), axG(2), 5,1, 4,2, 'S', 0, 'M', M, 'PT', 0.03, 'PB', 0.02);
    axAbsOrientCam =    subaxis(axG(1), axG(2), 5,3, 4,2, 'S', 0, 'M', M, 'PT', 0.02, 'PB', 0.02);

    axdRotLab =         subaxis(axG(1), axG(2), 9,1, 4,2, 'S', 0, 'M', M, 'PT', 0.03, 'PB', 0.02);
    axAbsOrientLab =    subaxis(axG(1), axG(2), 9,3, 4,2, 'S', 0, 'M', M, 'PT', 0.02, 'PB', 0.02);
    
    axIntXY =           subaxis(axG(1), axG(2), 5,5, 4,2, 'S', 0, 'M', M, 'PT', 0.02, 'PB', 0.03);
    axFWSideMove =      subaxis(axG(1), axG(2), 9,5, 4,2, 'S', 0, 'M', M, 'PT', 0.02, 'PB', 0.03);
    
    
    % Movie frame plot
    imshow(currFrame, [], 'Parent', axVid);
    axis off
    axVid.XTickLabel = [];
    
    % Movement map plot
    axes(ax2DPos);
    hold on
    nSegs = trialDuration;
    vectorLen = floor(size(mmXY, 1) / nSegs);
    cm = jet(nSegs);
    for iSeg = 1:nSegs
        if iSeg == 1
            pad = 1;
        else
            pad = 0;
        end
        currSegInd = (pad + ((iSeg-1)*vectorLen)):(iSeg * vectorLen);
        plot(ax2DPos, mmXY(currSegInd, 1), mmXY(currSegInd, 2), 'Color', cm(iSeg, :))
    end
    axis equal; %axis square
    lims = 1.1 * max(abs([ylim(ax2DPos), xlim(ax2DPos)]));
    if lims < 1
        lims = 1;
    end
    xlim([-2 7])
    ylim([-2 7])
    
    legend({'2D movement (rad)'}, 'FontSize', 8)
    x = mmXY(iFrame, 1);
    y = mmXY(iFrame, 2);
    [arrowVec(1), arrowVec(2)] = pol2cart(HD(iFrame)+(pi * 1.5), 0.5);
    %         [testVec(1), testVec(2)] = pol2cart(test(iFrame) + (pi * 1.5), 0.5);
    arrow([x - arrowVec(1)/2, y-arrowVec(2)/2], [x + arrowVec(1)/2, y + arrowVec(2)/2], 'FaceColor', 'green');
    %         arrow([x - testVec(1)/2, y-testVec(2)/2], [x + testVec(1)/2, y + testVec(2)/2], 'FaceColor', 'red');
    
    
    % Move direction plot
    axes(axMoveDirLab)
    plt = polarplot([moveDirLab(iFrame), moveDirLab(iFrame)], [0 1], '-', 'LineWidth', 5, 'Color', 'k');
    axMoveDirLab = plt.Parent;
    axMoveDirLab.FontSize = 10;
    axMoveDirLab.ThetaZeroLocation = 'bottom';
    axMoveDirLab.RTickLabel = [];
    axMoveDirLab.ThetaTickLabel = [];
    axMoveDirLab.RTick = [];
    axes(axMoveDirLab);
    hold on
    polarplot(linspace(0, pi * 2, 100), ones(100,1), 'linewidth', 0.5, 'Color', 'k');
    title('MoveDir')
    
    % Head direction plot
    axes(axHD)
    plt = polarplot([HD(iFrame), HD(iFrame)], [0 1], '-', 'LineWidth', 5, 'Color', 'k');
    axHD = plt.Parent;
    axHD.FontSize = 10;
    axHD.ThetaZeroLocation = 'bottom';
    axHD.RTickLabel = [];
    axHD.RTick = [];
    axHD.ThetaTickLabel = [];
    axes(axHD);
    hold on
    polarplot(linspace(0, pi * 2, 100), ones(100,1), 'linewidth', 0.5, 'Color', 'k');
    title('int HD')
    
    % dRotCam plot
    axes(axdRotCam);
    hold on
    plot(dRotCam);
    axdRotCam.XTick = xTickFR(1:5:end);
    axdRotCam.XTickLabel = []%xTickLabels(1:5:end);
    legend({'1', '2', '3'}, 'FontSize', 8, 'Location', 'NorthEast')
    plot(iFrame, dRotCam(iFrame, 1), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    plot(iFrame, dRotCam(iFrame, 2), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    plot(iFrame, dRotCam(iFrame, 3), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    title('dRotCam')
    
    % absOrientCam plot
    axes(axAbsOrientCam);
    hold on
    plot(absOrientCam);
    axAbsOrientCam.XTick = []%xTickFR(1:5:end);
    axAbsOrientCam.XTickLabel = xTickLabels(1:5:end);
    legend({'1', '2', '3'}, 'FontSize', 8, 'Location', 'NorthEast')
    plot(iFrame, absOrientCam(iFrame, 1), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    plot(iFrame, absOrientCam(iFrame, 2), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    plot(iFrame, absOrientCam(iFrame, 3), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    title('absOrientCam')
    
    % dRotLab plot
    axes(axdRotLab);
    hold on
    plot(dRotLab);
    axdRotLab.XTick = xTickFR(1:5:end);
    axdRotLab.XTickLabel = []%xTickLabels(1:5:end);
    legend({'1', '2', '3'}, 'FontSize', 8, 'Location', 'NorthEast')
    plot(iFrame, dRotLab(iFrame, 1), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    plot(iFrame, dRotLab(iFrame, 2), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    plot(iFrame, dRotLab(iFrame, 3), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    title('dRotLab')
    
    % absOrientLab plot
    axes(axAbsOrientLab);
    hold on
    plot(absOrientLab);
    axAbsOrientLab.XTick = xTickFR(1:5:end);
    axAbsOrientLab.XTickLabel = []%xTickLabels(1:5:end);
    legend({'1', '2', '3'}, 'FontSize', 8, 'Location', 'NorthEast')
    plot(iFrame, absOrientLab(iFrame, 1), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    plot(iFrame, absOrientLab(iFrame, 2), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    plot(iFrame, absOrientLab(iFrame, 3), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    title('absOrientLab')
    
   % XY displacement plot
    axes(axIntXY);
    hold on
    plot(mmXY);
    axIntXY.XTick = xTickFR(1:5:end);
    axIntXY.XTickLabel = xTickLabels(1:5:end);
    legend({'1', '2'}, 'FontSize', 8, 'Location', 'NorthEast')
    plot(iFrame, mmXY(iFrame, 1), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    plot(iFrame, mmXY(iFrame, 2), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    title('int XY displacement')
      
   % FW/side move plot
    axes(axFWSideMove);
    hold on
    plot(intFWSideMove);
    axFWSideMove.XTick = xTickFR(1:5:end);
    axFWSideMove.XTickLabel = xTickLabels(1:5:end);
    legend({'1', '2'}, 'FontSize', 8, 'Location', 'NorthWest')
    plot(iFrame, intFWSideMove(iFrame, 1), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    plot(iFrame, intFWSideMove(iFrame, 2), 'p', 'Color', 'k', 'MarkerSize', markerSize, 'MarkerFaceColor', 'k')
    title('FW/Side move')  
    axis on
    
    
    % Write frame to video(s)
    writeFrame = getframe(h);
    writeVideo(myVidWriter, writeFrame);
end%iFrame

close(myVidWriter)
