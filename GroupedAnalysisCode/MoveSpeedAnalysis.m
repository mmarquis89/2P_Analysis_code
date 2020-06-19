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
        flData = []; rawFlData = [];
        allVolTimes = [];
        sData = obj.sourceDataTable.subset;
        speedMat = generate_speedMat(obj);
        flMat = generate_flMat(obj);
        rawFlMat = generate_flMat(obj, 'rawFl');
        
        if ~isempty(obj.params.skipTrials)
            sData(obj.params.skipTrials, :) = [];
            speedMat(:, obj.params.skipTrials) = [];
            flMat(:, obj.params.skipTrials) = [];
            rawFlMat(:, obj.params.skipTrials) = [];
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
            currRawFlData = rawFlMat(:, iTrial);
            
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
                shutoffVols = tData.pmtShutoffVols{:};
                if max(shutoffVols) > numel(currFlData)
                   disp(['PMT shutoff vols: ', num2str(shutoffVols)])
                   disp(['Trial ', num2str(iTrial), ' nVolumes = ' num2str(numel(currFlData))]); 
                end
                shutoffVols = shutoffVols(shutoffVols <= numel(currFlData));
                speedDataVols(shutoffVols) = nan;
                currFlData(shutoffVols) = nan;
            end
            
            % Exclude event volumes as neccessary
            for iType = 1:numel(obj.eventNames)
                currEvent = obj.eventNames{iType};
                if obj.params.(currEvent)
                    currEventData = innerjoin(tData(:, {'expID', 'trialNum'}), ...
                        obj.eventObjects.(currEvent).eventData);
                    
                    for iEvent = 1:size(currEventData, 1)
                        if strcmp(currEvent, 'locomotion')
                            % Get rid of volumes on either side of the onset or offset of an event
                            currOnsetTime = currEventData(iEvent, :).onsetTime;
                            currOffsetTime = currEventData(iEvent, :).offsetTime;
                            excludeDur = obj.params.(currEvent);
                            onsetExcludeVols = (volTimes > (currOnsetTime - excludeDur)) & ...
                                    volTimes < (currOnsetTime + excludeDur);
                            offsetExcludeVols = (volTimes > (currOffsetTime - excludeDur)) & ...
                                    volTimes < (currOffsetTime + excludeDur);
                            speedDataVols(find(onsetExcludeVols | offsetExcludeVols)) = nan;
                            currFlData(find(onsetExcludeVols | offsetExcludeVols)) = nan;
                            
                            
                        else
                            % Get rid of volumes during or just after an event (i.e. sensory stim)
                            excludeVols = (volTimes > currEventData(iEvent, :).onsetTime & ...
                                volTimes < (currEventData(iEvent, :).offsetTime + ...
                                obj.params.(currEvent)));
                            speedDataVols(find(excludeVols)) = nan; %#ok<*FNDSB>
                            currFlData(find(excludeVols)) = nan;
                            
                            
                        end
                    end
                end                
            end
            
            % Calculate cross-correlation between Fl data and moveSpeed
            [r(iTrial, :), lags(iTrial, :)] = xcorr(currFlData(~isnan(currFlData + speedDataVols)), ...
                    speedDataVols(~isnan(currFlData + speedDataVols)), 6); 
            
            % Apply a lag to the speed data if necessary
            speedDataVols = speedDataVols(1:end - obj.params.lagVols);
            currFlData = currFlData((obj.params.lagVols + 1):end);
            currRawFlData = currRawFlData((obj.params.lagVols + 1):end);
            volTimes = volTimes(1:end - obj.params.lagVols);
            % Add to plotting vectors
            speedData = [speedData; speedDataVols];
            flData = [flData; currFlData];
            rawFlData = [rawFlData; currRawFlData];
            allVolTimes = [allVolTimes; volTimes' + (tData.trialDuration * (iTrial - 1))];
        end
        
        % Drop volumes where either fl or speed data is NaN
        nanVols = isnan(flData) | isnan(speedData) | flData == 0 | speedData == 0;
        speedData(nanVols) = [];
        flData(nanVols) = [];
        rawFlData(nanVols) = [];
%         allVolTimes(nanVols) = [];
        
        % Save processed data
        output = struct();
        output.speedData = speedData;
        output.flData = flData;
        output.rawFlData = rawFlData;
        output.volTimes = allVolTimes;
        output.r = r;
        output.lags = lags;
        output.binnedData = [];
        output.paramSnapshot = obj.params;
        obj.analysisOutput = output;
        
    end
    
    % Calculate mean fl value in equal-sized bins of analyzed moveSpeed data
    function obj = generate_binned_flData(obj, flType, speedNorm)
        aData = obj.analysisOutput;
        if ~isempty(obj.params.slidingWinSize)
            slidingWinSize = obj.params.slidingWinSize;
        else
            slidingWinSize = floor(numel(aData.rawFlData) / obj.params.nAnalysisBins);
        end
        if nargin < 2
            flType = obj.analysisOutput.paramSnapshot.flType;
        end
        if strcmpi(flType, 'slidingBaseline')
            base = mov_percentile(aData.rawFlData, slidingWinSize, 0.05)';
            flData = (aData.rawFlData - base) ./ base;
        elseif strcmpi(flType, 'normalized')
            flData = mov_norm(aData.rawFlData, slidingWinSize)';
        else
            flData = aData.flData;
        end
        if nargin == 3 && speedNorm == 1
            speedData = mov_norm(aData.speedData, slidingWinSize)';  
        else
            speedData = aData.speedData;
            speedNorm = 0;
        end
        
        
        [speedBinMidpoints, flBinMeans, SEM] = obj.bin_data(speedData, flData, ...
                obj.params.nAnalysisBins, obj.params.plotting_minSpeed);

        obj.analysisOutput.binnedData.speedBinMidpoints = speedBinMidpoints;
        obj.analysisOutput.binnedData.flBinMeans = flBinMeans;
        obj.analysisOutput.binnedData.flBinSEM = SEM;
        obj.analysisOutput.binnedData.flType = flType;
        obj.analysisOutput.binnedData.speedNorm = speedNorm;
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
                '-o', 'color', 'k', 'linewidth', 2);
        hold on;
% yyaxis right
% a = obj.analysisOutput.corrTestFl(:);
% b = obj.analysisOutput.corrTestSpeed(:);
% c = isnan(a) | isnan(b);
% [r, lags] = xcorr(a(~c), b(~c), 6);
% plot(ax, lags, r, '-o', 'color', 'r', 'linewidth', 2);
% legend({'trial-averaged', 'concatenated'}, 'autoupdate', 'off', 'location', 'nw')
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
        volTimes = obj.analysisOutput.volTimes;
        medVolDur = median(diff(volTimes), 'omitnan');
        plotVolTimes = (1:numel(speedData)) * medVolDur;
        plot(plotVolTimes, speedData, 'color', 'b');
        ylabel('Speed (mm/sec)')
        xlabel('Time (sec)')
        hold on;
        plot(plotVolTimes, ones(size(plotVolTimes)), 'color', 'k')
        yyaxis right
        plot(plotVolTimes, flData, 'color', 'r');
        ylabel(obj.analysisOutput.paramSnapshot.flType);
        if nargin < 2 || isempty(nSamples)
            nSamples = numel(speedData);
        end
        xlim([0 nSamples]);
        legend({'moveSpeed', obj.analysisOutput.paramSnapshot.flType})
    end
    
    % ---------- Generate primary analysis plots ----------
    function ax = plot_2D_hist(obj, ax, flType, speedNorm)
        if nargin < 2
            f = figure(6); clf;
            f.Color = [1 1 1];
            ax = axes();
        end
        aData = obj.analysisOutput;
        if nargin < 3 || isempty(flType)
            flType = aData.paramSnapshot.flType;
        end
        if ~isempty(obj.params.slidingWinSize)
            slidingWinSize = obj.params.slidingWinSize;
        else
            slidingWinSize = floor(numel(aData.rawFlData) / obj.params.nAnalysisBins);
        end
        if strcmpi(flType, 'slidingBaseline')
            base = mov_percentile(aData.rawFlData, slidingWinSize, 0.05)';
            flData = (aData.rawFlData - base) ./ base;
        elseif strcmpi(flType, 'normalized')
            flData = mov_norm(aData.rawFlData, slidingWinSize)';
        else
            flData = aData.flData;
        end
        if nargin == 4 && speedNorm == 1
            speedData = mov_norm(aData.speedData, slidingWinSize)';
        else
            speedData = aData.speedData;
            speedNorm = 0;
        end
        speedData(speedData < obj.params.plotting_minSpeed) = nan;
        flData(speedData < obj.params.plotting_minSpeed) = nan;
        h = histogram2(ax, speedData, flData, obj.params.nHistBins, 'displaystyle', 'tile');
        xEdges = h.XBinEdges; 
        yEdges = h.YBinEdges;
        binCounts = h.BinCounts;
        if ~isempty(obj.params.maxHistBinCount)
            binCounts(binCounts > obj.params.maxHistBinCount) = obj.params.maxHistBinCount;
        end
        histogram2(ax, 'XBinEdges', xEdges, 'YBinEdges', yEdges, 'BinCounts', ...
                binCounts, 'displaystyle', 'tile');
        ax.XGrid = 'off';
        ax.YGrid = 'off';
        if speedNorm
            xlabel('Normalized moveSpeed');
        else
            xlabel('moveSpeed (mm/sec)')
        end
        if strcmpi(flType, 'normalized')
            ylabel('Normalized F');
        elseif strcmpi(flType, 'slidingbaseline')
            ylabel('Sliding dF/F')
        else
            ylabel(flType);
        end
    end
    function ax = plot_binned_fl(obj, ax)
        if nargin < 2
            f = figure(7); clf;
            f.Color = [1 1 1];
            ax = axes();
        end    
        binnedData = obj.analysisOutput.binnedData;
        errorbar(ax, binnedData.speedBinMidpoints, binnedData.flBinMeans, binnedData.flBinSEM, ...
                '-o', 'color', 'k', 'linewidth', 1);
        if binnedData.speedNorm
            xlabel('Normalized moveSpeed');
        else
            xlabel('moveSpeed (mm/sec)')
        end
        if strcmpi(binnedData.flType, 'normalized')
            ylabel('Mean normalized F');
        else
            ylabel(['Mean ', binnedData.flType]);
        end
    end
end%methods  

methods(Static)
    % Separate fluorescence data sorted by moveSpeed into equal-size bins and average across them
    function [binMidpoints, binMeans, binSEM] = bin_data(xData, yData, nBins, minX)
        if nargin < 4
           minX = -100; 
        end
        plotDataTable = table(xData, yData, 'variableNames', {'speed', 'fl'});
        plotDataSort = sortrows(plotDataTable, 'speed');
        plotDataSort = plotDataSort(plotDataSort.speed > minX, :);
        binSize = round(size(plotDataSort, 1) / nBins);
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
        binMidpoints = speedBinMidpoints;
        binMeans = flBinMeans;
        binSEM = SEM;
    end
end


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
    paramStruct.maxHistBinCount = [];
    paramStruct.nAnalysisBins = 100; 
    paramStruct.plotting_minSpeed = 0;
    paramStruct.skipTrials = [];
    paramStruct.nAnalysisBins = 30;
    paramStruct.flDataType = 'rawFl';
    paramStruct.slidingWinDur = 60;
    paramStruct.slidingWinSize = [];
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
    outputMat = cell2padded_mat(currSubset.fwSpeed); % --> [frame, trial]
    outputMat = outputMat * 25 * 4.5; % Convert from rad/frame to mm/sec
    outputMat = repeat_smooth(outputMat, p.nSmoothReps, 'dim', 1, 'smWin', p.smWinFrames);
    outputMat(outputMat > p.max_moveSpeed) = nan;
end

% Generate a matrix containing smoothed and padded fluorescence data for all trials in subset
function outputMat = generate_flMat(obj, flType)
    if nargin < 2
       flType = obj.params.flType; 
    end
    currSubset = obj.sourceDataTable.subset;
    outputMat = cell2padded_mat(currSubset.rawFl); % --> [vol, trial]
    outputMat = smoothdata(outputMat, 1, 'gaussian', obj.params.smWinVols, 'omitnan');
    if strcmp(flType, 'trialDff')
        baseF = repmat(obj.sourceDataTable.subset.trialBaseline', size(outputMat, 1), 1);
        outputMat = (outputMat - baseF) ./ baseF;
    elseif strcmp(flType, 'expDff')
        baseF = repmat(obj.sourceDataTable.subset.expBaseline', size(outputMat, 1), 1);
        outputMat = (outputMat - baseF) ./ baseF;
    end
end



