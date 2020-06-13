classdef MoveSpeedAnalysis
% ==================================================================================================   
%    
% Properties:
%
%       sourceData
%       eventObjects
%       filterDefs
%       params
%       
% Methods:
%
% 
%
% ==================================================================================================    
    
properties
    sourceDataTable % DataTable with all data except events
    eventObjects
    filterDefs
    params
    analysisOutput
end

properties (Dependent)
    eventNames
end

methods
    
    % Constructor
    function obj = MoveSpeedAnalysis(expList)
        
        % Load metadata, FicTrac data, ROI data for all experiments
        parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';
        allExpData = [];
        for iExp = 1:size(expList, 1)
            
            currExpID = expList.expID{iExp};
            disp(currExpID);
            
            % Generate full paths to data files
            roiFileName = [currExpID, '_roiData.mat'];
            roiFile = fullfile(parentDir, roiFileName);
            ftFileName = [currExpID, '_ficTracData.mat'];
            ftFile = fullfile(parentDir, ftFileName);
            expMdFile = fullfile(parentDir, [currExpID, '_expMetadata.csv']);
            trialMdFile = fullfile(parentDir, [currExpID, '_trialMetadata.mat']);
            
            if exist(roiFile, 'file') && exist(ftFile, 'file') && exist(expMdFile, 'file')
                load(roiFile, 'roiData');
                load(ftFile, 'ftData');
                expMd = readtable(expMdFile, 'delimiter', ',');
                load(trialMdFile, 'trialMetadata');
                trialMd = trialMetadata;
                allExpData = [allExpData; inner_join(expMd, trialMd, roiData, ftData)];
            else
                disp(['Skipping ', currExpID, ' due to one or more missing files']);
            end
        end%iExp
        
        % Create source DataTable
        obj.sourceDataTable = DataTable(allExpData);
        obj.filterDefs = struct();
        obj.filterDefs.roiName = [];
        obj.filterDefs.expID = [];
        
        % Load all event data objects
        obj.eventObjects = load_event_data(expList, parentDir);
        
        % Set initial analysis parameters
        obj.params = initialize_analysis_params(obj);
    end
    
    %  ---------- Dependent property methods ----------
    function eventNames = get.eventNames(obj)
        eventNames = fieldnames(obj.eventObjects)';
    end
    
    % Initialize filters on the source DataTable using obj.filterDefs
    function obj = init_filters(obj) 
        obj.sourceDataTable = obj.sourceDataTable.initialize_filters(obj.filterDefs);
    end
    
    % Run the analysis processing steps on the data from the current subset
    function obj = analyze(obj)
        
        speedData = [];
        flData = [];
        
        sData = obj.sourceDataTable.subset;
        speedMat = generate_speedMat(obj);
        flMat = generate_flMat(obj);
        
        if ~isempty(obj.params.skipTrials)
            sData(obj.params.skipTrials, :) = [];
            speedMat(:, obj.params.skipTrials) = [];
            flMat(:, obj.params.skipTrials) = [];
        end
        
        r = []; lags = [];
        for iTrial = 1:size(sData, 1)
            if ~mod(iTrial, 10)
                disp([num2str(iTrial), ' of ', num2str(size(sData, 1))])
            end
            
            % Get data for current trial
            tData = sData(iTrial, :);
            speedDataFrames = speedMat(:, iTrial);
            currFlData = flMat(:, iTrial);
            
            % Calculate volume times
            volTimes = calc_volTimes(tData.nVolumes, tData.volumeRate, tData.trialDuration, ...
                tData.originalTrialCount);
            
            % Downsample moveSpeed data to match volTimes
            speedDataVols = [];
            for iVol = 1:numel(volTimes)
                speedDataVols(iVol) = speedDataFrames(argmin(abs(tData.frameTimes{:} - ...
                        volTimes(iVol))));
            end
            speedDataVols = speedDataVols';
            speedDataVols = [speedDataVols; nan(numel(currFlData) - numel(speedDataVols), 1)];
            
            % Exclude PMT shutoff vols
            if ~isempty(tData.pmtShutoffVols{:}) && ~any(isnan(tData.pmtShutoffVols{:}))
                speedDataVols(tData.pmtShutoffVols{:}) = nan;
                currFlData(tData.pmtShutoffVols{:}) = nan;
            end
            
            % Exclude event volumes as neccessary
            for iType = 1:numel(obj.eventNames)
                currEvent = obj.eventNames{iType};
                if obj.params.(currEvent)
                    currEventData = innerjoin(tData(:, {'expID', 'trialNum'}), ...
                        obj.eventObjects.(currEvent).eventData);
                    for iEvent = 1:size(currEventData, 1)
                        excludeVols = (volTimes > currEventData(iEvent, :).onsetTime & volTimes < ...
                                (currEventData(iEvent, :).offsetTime + obj.params.(currEvent)));
                        speedDataVols(find(excludeVols)) = nan; %#ok<*FNDSB>
                        currFlData(find(excludeVols)) = nan;
                    end
                end                
            end
            
            % Calculate cross-correlation between Fl data and moveSpeed
            [r(iTrial, :), lags(iTrial, :)] = xcorr(currFlData(~isnan(currFlData + speedDataVols)), ...
                    speedDataVols(~isnan(currFlData + speedDataVols)), 6); 
            
            % Apply a lag to the speed data if necessary
            speedDataVols = speedDataVols(1:end - obj.params.lagVols);
            currFlData = currFlData((obj.params.lagVols + 1):end);
            
            % Add to plotting vectors
            speedData = [speedData; speedDataVols];
            flData = [flData; currFlData];
        end
        
        % Drop volumes where either fl or speed data is NaN
        nanVols = isnan(flData) | isnan(speedData);
        speedData(nanVols) = [];
        flData(nanVols) = [];
        
        % Save processed data
        output = struct();
        output.speedData = speedData;
        output.flData = flData;
        output.r = r;
        output.lags = lags;
        output.speedBinMidpoints = [];
        output.flBinMeans = [];
        output.flBinSEM = [];
        output.paramSnapshot = obj.params;
        obj.analysisOutput = output;
        
    end
    
    % Calculate mean fl value in equal-sized bins of analyzed moveSpeed data
    function obj = generate_binned_flData(obj)
        aData = obj.analysisOutput;
        plotDataTable = table(aData.speedData, aData.flData, 'variableNames', {'speed', 'fl'});
        plotDataSort = sortrows(plotDataTable, 'speed');
        plotDataSort = plotDataSort(plotDataSort.speed > obj.params.plotting_minSpeed, :);
        binSize = round(size(plotDataSort, 1) / obj.params.nAnalysisBins);
        binEdges = 1:binSize:(size(plotDataSort, 1) - 1);
        flBinMeans = []; SEM = []; speedBinMidpoints = [];
        for iBin = 1:numel(binEdges)
            if iBin < numel(binEdges)
                binInds = binEdges(iBin):(binEdges(iBin + 1) - 1);
            else
                binInds = binEdges(iBin):size(plotDataSort, 1);
            end
            flBinMeans(iBin) = mean(plotDataSort{binInds, 'fl'}, 'omitnan');
            SEM(iBin) = std_err(plotDataSort{binInds, 'fl'});
            speedBinMidpoints(iBin) = plotDataSort{binInds(round(numel(binInds) / 2)), 'speed'};
        end
        obj.analysisOutput.speedBinMidpoints = speedBinMidpoints;
        obj.analysisOutput.flBinMeans = flBinMeans;
        obj.analysisOutput.flBinSEM = SEM;
    end
    
    % ---------- Generate summary plots (pre-analysis) ----------
    function ax = plot_speedMat(obj, ax)
        if nargin < 2
            f = figure(1); clf;
            f.Color = [1 1 1];
            ax = axes();
        end
        speedMat = generate_speedMat(obj); % --> [frame, trial]
        if ~isempty(obj.params.skipTrials)
            speedMat(:, obj.params.skipTrials) = [];
        end
        imagesc(ax, speedMat');
        colorbar()
        title('MoveSpeed (mm/sec)')
        xlabel('Frame')
        ylabel('Trial')
    end
    function ax = plot_sortedSpeed(obj, ax)
        if nargin < 2
            f = figure(2); clf;
            f.Color = [1 1 1];
            ax = axes();
        end
        speedMat = generate_speedMat(obj); % --> [frame, trial]
        if ~isempty(obj.params.skipTrials)
            speedMat(:, obj.params.skipTrials) = [];
        end
        plot(ax, sort(speedMat(:)), 'color', 'k', 'linewidth', 2)
        hold on;
        xL = [0, numel(speedMat(~isnan(speedMat)))];
        yy = [1 1] * obj.params.plotting_minSpeed; 
        plot(xL, yy, 'color', 'r', 'linewidth', 2)
        xlim(xL)
        ylabel('Sorted moveSpeed (mm/Sec)')
    end  
    function ax = plot_flMat(obj, ax)
        if nargin < 3
            f = figure(3); clf;
            f.Color = [1 1 1];
            ax = axes();
        end
        flType = obj.params.flType;
        flMat = generate_flMat(obj);
        if ~isempty(obj.params.skipTrials)
            flMat(:, obj.params.skipTrials) = [];
        end
        if strcmp(flType, 'trialDff')
            titleStr = 'Trial dF/F';
        elseif strcmp(flType, 'expDff')
            titleStr = 'Full exp dF/F';
        else
            titleStr = 'Raw Fl';
        end
        imagesc(ax, flMat');
        title(titleStr);
        colorbar()
        xlabel('Volume')
        ylabel('Trial')
    end
    
    % ---------- Generate post-analysis summary plots ----------
    function ax = xcorr_plot(obj, ax)
        if nargin < 2
            f = figure(4); clf;
            f.Color = [1 1 1];
            ax = axes();
        end
        plot(ax, obj.analysisOutput.lags(1, :), mean(obj.analysisOutput.r, 1, 'omitnan'), ...
                '-o', 'color', 'k', 'linewidth', 3);
        hold on;
        yL = ylim();
        plot([0 0], yL, '--', 'color', 'b')
        ax.XAxisLocation = 'top';
        xlabel('Lag (volumes)')
        ylabel ('r value')
    end
    function ax = norm_data_overlay(obj, nSamples, ax)
        if nargin < 3
            f = figure(5); clf;
            f.Color = [1 1 1];
            ax = axes();
        end
        speedData = obj.analysisOutput.speedData;
        flData = obj.analysisOutput.flData;
        plot(speedData, 'color', 'b');
        ylabel('Speed (mm/sec)')
        hold on;
        yyaxis right
        plot(flData, 'color', 'r');
        ylabel(obj.analysisOutput.paramSnapshot.flType);
        if nargin < 2 || isempty(nSamples)
            nSamples = numel(speedData);
        end
        xlim([0 nSamples]);
        legend({'moveSpeed', obj.analysisOutput.paramSnapshot.flType})
    end
    
    % ---------- Generate primary analysis plots ----------
    function ax = plot_2D_hist(obj, ax)
        if nargin < 2
            f = figure(6); clf;
            f.Color = [1 1 1];
            ax = axes();
        end
        aData = obj.analysisOutput;
        fl = aData.flData;
        speed = aData.speedData;
        speed(speed < obj.params.plotting_minSpeed) = nan;
        fl(speed < obj.params.plotting_minSpeed) = nan;
        histogram2(ax, speed, fl, obj.params.nHistBins, 'displaystyle', 'tile');
        ax.XGrid = 'off';
        ax.YGrid = 'off';
        xlabel('moveSpeed (mm/sec)'); 
        ylabel(aData.paramSnapshot.flType); 
    end
    function ax = plot_binned_fl(obj, ax)
        if nargin < 2
            f = figure(7); clf;
            f.Color = [1 1 1];
            ax = axes();
        end
        aData = obj.analysisOutput;
        errorbar(ax, aData.speedBinMidpoints, aData.flBinMeans, aData.flBinSEM, '-o', 'color', 'k');
        xlabel('moveSpeed (mm/sec)')
        ylabel(['mean binned ', aData.paramSnapshot.flType]);
    end
    
    
end%methods

end%class

% ==================================================================================================
% Local functions
% ==================================================================================================
function paramStruct = initialize_analysis_params(obj)
    paramStruct = struct();
    paramStruct.max_moveSpeed = 100;
    paramStruct.smWinVols = 1;
    paramStruct.smWinFrames = 1;
    paramStruct.nSmoothRepsFrames = 15;
    paramStruct.lagVols = 2;
    paramStruct.nHistBins = 50;
    paramStruct.nAnalysisBins = 100; 
    paramStruct.plotting_minSpeed = 0;
    paramStruct.skipTrials = [];
    paramStruct.flDataType = 'rawFl';
    paramStruct.excludeEvents = struct();
    for iType = 1:numel(obj.eventNames)
        paramStruct.(obj.eventNames{iType}) = 0;
    end
end

% Takes a cell array of variable-length column vectors and converts them to a matrix by padding all
% vectors with NaN to match the size of the longest one
function outputMat = cell2padded_mat(inputCell)
    maxLen = max(cellfun(@numel, inputCell));
    outputMat = [];
    for i = 1:numel(inputCell)
        if numel(inputCell{i}) == maxLen
            outputMat(:, i) = inputCell{i};
        else
            padVec = nan(maxLen - numel(inputCell{i}), 1);
            outputMat(:, i) = [inputCell{i}; padVec];
        end        
    end
end

% Generate a matrix containing smoothed, padded, and capped moveSpeed data for all trials in subset
function outputMat = generate_speedMat(obj)
    currSubset = obj.sourceDataTable.subset;
    p = obj.params;
    outputMat = cell2padded_mat(currSubset.moveSpeed); % --> [frame, trial]
    outputMat = outputMat * 25 * 4.5; % Convert from rad/frame to mm/sec
    outputMat = repeat_smooth(outputMat, p.nSmoothReps, 'dim', 1, 'smWin', p.smWinFrames);
    outputMat(outputMat > p.max_moveSpeed) = p.max_moveSpeed;
end

% Generate a matrix containing smoothed and padded fluorescence data for all trials in subset
function outputMat = generate_flMat(obj)
    currSubset = obj.sourceDataTable.subset;
    outputMat = cell2padded_mat(currSubset.rawFl); % --> [vol, trial]
    outputMat = smoothdata(outputMat, 1, 'gaussian', obj.params.smWinVols, 'omitnan');
    if strcmp(obj.params.flType, 'trialDff')
        baseF = repmat(obj.sourceDataTable.subset.trialBaseline', size(outputMat, 1), 1);
        outputMat = (outputMat - baseF) ./ baseF;
    elseif strcmp(obj.params.flType, 'expDff')
        baseF = repmat(obj.sourceDataTable.subset.expBaseline', size(outputMat, 1), 1);
        outputMat = (outputMat - baseF) ./ baseF;
    end
end









