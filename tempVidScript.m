%% CREATE FICTRAC + BEHAVIOR VIDS
parentDir = 'B:\Dropbox (HMS)\2P Data\Behavior Vids\2018_04_05_exp_2\_Movies'

vidFiles = dir(fullfile(parentDir, ['sid*tid*.mp4']));

for iTrial = 1:size(selectData, 3)
    
    vidName = vidFiles(iTrial).name;
    disp(vidName)
    vidFile = fullfile(parentDir, vidName);
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
        axDisp = subaxis(3,6,[13:15], 'S', 0, 'M', M, 'PB', 0.05, 'PR', 0.01, 'PL', 0.008);
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
        lims = max(abs([ylim(axMove), xlim(axMove)]));
        xlim([-lims lims])
        ylim([-lims lims])
        axis equal; %axis square
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
        axes(axDisp);
        hold on
        plot(mmData(:, [3 4]));
        axDisp.XTick = xTick;
        axDisp.XTickLabel = xTickLabels;
        legend({'X disp (mm)', 'Y disp (mm)'}, 'FontSize', 11, 'Location', 'NorthWest')
        xlabel('Time (sec)')
        plot(iFrame, mmData(iFrame, 3), 'p', 'Color', 'k', 'MarkerSize', 12, 'MarkerFaceColor', 'k')
        plot(iFrame, mmData(iFrame, 4), 'p', 'Color', 'k', 'MarkerSize', 12, 'MarkerFaceColor', 'k')
        
        % Velocity plot
        axes(axVel)
        plot(velData(:,5));
        axVel.XTick = xTick;
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

%% CREATE VID WITH TRIAL-BY-TRIAL FICTRAC PLOTS

% Create VidWriter
myVidWriter = VideoWriter(fullfile(parentDir, 'FicTrac_Trial_Plots.avi'));
myVidWriter.FrameRate = 1;
open(myVidWriter)

for iTrial = 1:size(selectData, 3)
    
    disp(num2str(iTrial))
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
    
    % Create figure
    h = figure(100); clf
    screenSize = get( groot, 'Screensize' );
    h.Color = [1 1 1];
    h.OuterPosition = [50 100 (screenSize(3) - 100) screenSize(4) - 150];
    
    % Create axes
    M = 0.01;
    P = 0.00;
    axMove = subaxis(3,6,[1:3 7:9 13:15], 'S', 0, 'M', M, 'PB', 0.03);
    axHD = subaxis(3,6,[4:6], 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.00);
    axDisp = subaxis(3,6,[10:12], 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.00);
    axVel = subaxis(3,6,[16:18], 'S', 0, 'M', M, 'PB', 0.05, 'PL', 0.00);
    
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
        plot(axMove, smooth(mmData(currSegInd, 3), 3), smooth(mmData(currSegInd, 4), 3), 'Color', cm(iSeg, :))
    end
    axis equal;
    lims = 1.1 * (max(abs([ylim(axMove), xlim(axMove)])));
    if lims < 1
        lims = 1;
    end
    xlim([-lims lims])
    ylim([-lims lims])
    legend({'2D movement (mm)'}, 'FontSize', 11)
    x = mmData(iFrame, 3);
    y = mmData(iFrame, 4);
    
    % Head direction plot
    axes(axHD)
    plot(HD);
    legend({'Heading (rad)'}, 'FontSize', 11, 'Location', 'NorthWest');
    axHD.XTick = xTick;
    axHD.XTickLabel = xTickLabels;
    xlabel('Time (sec)')
    
    % Displacement plot
    axes(axDisp);
    hold on
    plot(mmData(:, [3 4]));
    axDisp.XTick = xTick;
    axDisp.XTickLabel = xTickLabels;
    legend({'X disp (mm)', 'Y disp (mm)'}, 'FontSize', 11, 'Location', 'NorthWest')
    xlabel('Time (sec)')
    
    % Velocity plot
    axes(axVel)
    plot(velData(:,5));
    axVel.XTick = xTick;
    axVel.XTickLabel = xTickLabels;
    legend({'Speed (mm/sec)'}, 'FontSize', 11)
    xlabel('Time (sec)')
    hold on

    % Write frame to video(s)
    writeFrame = getframe(h);
    writeVideo(myVidWriter, writeFrame);
    
end

close(myVidWriter)
