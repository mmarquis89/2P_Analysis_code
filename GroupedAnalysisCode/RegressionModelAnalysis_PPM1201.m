classdef RegressionModelAnalysis_PPM1201
% ==================================================================================================
%
% Properties:
%       sourceData (table)
%       sourceDataParams (scalar struct)
%       modelData (table)
%       modelParams (scalar struct)
%
% Methods:
%       initialize_models(modelParams)
%       optimize_odor_integration_windows()
%       optimize
%       ax = plot_odor_filter(expID, axesHandle)       
%       plot_mean_moveSpeed(expID, axesHandle)
%       
%       (Static) plot_coeffs(mdl, axesHandle);
% ==================================================================================================

properties
    sourceData
    sourceDataParams
    modelData
    modelParams
end

methods
    
    % Constructor
    function obj = RegressionModelAnalysis_PPM1201(expInfoTbl, sourceDataParams)
        
        % expInfoTable = table with colums: [expID][skipTrials][skipVols]
        
        % sourceDataParams = struct with fields: maxSpeed, smWinVols, smWinFrames, smReps, 
        %                                        ftLagVols, speedType, odorRespFilterDur
        
        obj.sourceDataParams = sourceDataParams;
        
        % Unpack sourceDataParams
        maxSpeed = sourceDataParams.maxSpeed;
        smWinVols = sourceDataParams.smWinVols;
        smWinFrames = sourceDataParams.smWinFrames;
        smReps = sourceDataParams.smReps;
        ftLagVols = sourceDataParams.ftLagVols;
        roiName = sourceDataParams.roiName;
        obj.modelData = [];
        
        % -------  Load source data for all experiments -------
        parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
        dataDir = fullfile(parentDir, 'all_experiments');
        FRAME_RATE = 25;
        
        % Load AlignEvent object
        disp('Loading new AlignEvent object...');
        load(fullfile(parentDir, 'Saved_AlignEvent_objects', '20200609_AlignEventObj_odor.mat'), ...
                'alignObj');
        disp('Loading complete')
        
        % Load event DataTable
        disp('Loading event DataTable...')
        fileName = ['odor_pre_', num2str(sourceDataParams.odorRespFilterDur(1)), '_post_', ...
                num2str(sourceDataParams.odorRespFilterDur(2)), '.mat'];
        load(fullfile(parentDir, 'Saved_DataTable_objects', fileName), 'dt')
        disp('Loading complete')
        
        % Pre-process data for each experiment
        obj.sourceData = [];
        for iExp = 1:numel(expInfoTbl.expID)
            
            currExpID = expInfoTbl.expID{iExp};
            
            % -------  Load source data -------
            [expMd, trialMd] = load_metadata({currExpID}, dataDir);
            roiData = load_roi_data({currExpID}, dataDir);
            roiNames = regexp(roiData.roiName, [roiName, '-.'], 'match', 'once'); % Just one hemi
            roiNames(cellfun(@isempty, roiNames)) = [];   
            roiData = roiData(strcmp(roiData.roiName, roiNames{1}), :);
            ft = load_ft_data({currExpID}, dataDir);
            
            % Load odor event data
            eventData = load_event_data({currExpID}, dataDir);
            odorEventData = eventData.odor;
            
            % Process data for each trial (assuming no time between trials)
            trialNums = trialMd.trialNum;
            if ~isempty(expInfoTbl.skipTrials{iExp})
                trialNums(ismember(trialNums, expInfoTbl.skipTrials{iExp})) = [];
            end
            allFwSpeed = []; allYawSpeed = []; allOdorVols = []; allFl = []; allVolTimes = [];
            for iTrial = 1:numel(trialNums)
                currTrialNum = trialNums(iTrial);
                md = innerjoin(expMd, trialMd(currTrialNum, :));
                volTimes = calc_volTimes(md.nVolumes, md.volumeRate, md.trialDuration, ...
                        md.originalTrialCount);
                
                % Smooth imaging rawFl data for current trial
                currRoiData = roiData(currTrialNum, :);
                currFl = smoothdata(currRoiData.rawFl{:}, 'gaussian', smWinVols);
                if ~isnan(md.pmtShutoffVols{:})
                    currFl(md.pmtShutoffVols{:}) = nan;
                end
                
                % Smooth FicTrac data, then downsample to match volume rate
                currFt = ft(currTrialNum, :);
                if strcmp(sourceDataParams.speedType, 'moveSpeed')
                    currFwSpeed = repeat_smooth(currFt.moveSpeed{:} .* 4.5 .* FRAME_RATE, smReps, ...
                            'smWin', smWinFrames);
                else
                    currFwSpeed = repeat_smooth(currFt.fwSpeed{:} .* 4.5 .* FRAME_RATE, smReps, ...
                            'smWin', smWinFrames);
                end
                currFwSpeed(currFwSpeed > maxSpeed) = maxSpeed;
                currYawSpeed = abs(repeat_smooth(rad2deg(currFt.yawSpeed{:}) .* FRAME_RATE, ...
                        smReps, 'smWin', smWinFrames));
                currFwSpeedVols = []; currYawSpeedVols = [];
                for iVol = 1:numel(volTimes)
                    dsFrame = argmin(abs(currFt.frameTimes{:} - volTimes(iVol)));
                    currFwSpeedVols(iVol) = currFwSpeed(dsFrame);
                    currYawSpeedVols(iVol) = currYawSpeed(dsFrame);
                end
                
                % Apply a lag to the FicTrac data
                currFwSpeedVols = currFwSpeedVols([ones(1, ftLagVols), 1:end-ftLagVols]);
                currYawSpeedVols = currYawSpeedVols([ones(1, ftLagVols), 1:end-ftLagVols]);
                
                % Get odor command vector for current trial
                currTrialOdorVols = odorEventData.create_logical_array(md.nVolumes, volTimes, ...
                        table({currExpID}, currTrialNum, 'variableNames', {'expID', 'trialNum'}));
                
                % Append variables to whole-experiment vectors
                allFl = [allFl; currFl];
                allFwSpeed = [allFwSpeed; currFwSpeedVols'];
                allYawSpeed = [allYawSpeed; currYawSpeedVols'];
                allOdorVols = [allOdorVols; currTrialOdorVols];
                if iTrial == 1
                    allVolTimes = volTimes;
                else
                    allVolTimes = [allVolTimes, volTimes + allVolTimes(end)];
                end
            end%iTrial
            
            % Calculate the odor response filter from trial-averaged data
            disp('Generating odor response vector...');
            filterDefs = alignObj.create_filterDefs;
            filterDefs.expID = '20180329-2';
            filterDefs.roiName = 'TypeF-R';
            filterDefs.odor = -1;
            filterDefs.locomotion = 0;
            filterDefs.flailing = 0;
            filterDefs.grooming = 0;
            dt = dt.initialize_filters(filterDefs);

            meanTrace = smoothdata(mean(cell2mat(dt.subset.eventFl'), 2, 'omitnan'), ...
                    'gaussian', smWinVols);
            odorRespFilter = meanTrace(dt.subset.volTimes{1} > 0); % Trim pre-onset volumes
            odorRespFilter = odorRespFilter - odorRespFilter(1);   % Offset to zero first volume
            
            % Also calculate a trial-averaged version of the move speed while we're at it
            filterDefs.expID = currExpID;
            filterDefs.roiName = 'TypeF';
            filterDefs.locomotion = [];
            dt2 = dt.initialize_filters(filterDefs);
            moveSpeedArr = cell2mat(dt2.subset.moveSpeed') .* 4.5 .* FRAME_RATE;
            meanSpeed = repeat_smooth(mean(moveSpeedArr, 2, 'omitnan'), smReps, 'smWin', ...
                    smWinFrames);
            meanSpeed = meanSpeed(dt2.subset.frameTimes{1} > 0);
            
            % Extract onset volumes from odor command vector
            odorOnsetVols = (regexp(regexprep(num2str(allOdorVols'), ' ', ''), '01')) + 1;
            odorOnsetVector = zeros(size(allOdorVols));
            odorOnsetVector(odorOnsetVols) = 1;
            
            % Convolve with odor response filter to get predicted fl responses
            odorRespVector = conv(odorOnsetVector, odorRespFilter);
            odorRespVector = odorRespVector(1:numel(odorOnsetVector));
            
            % Add current experiment to sourceData table
            newRow = expInfoTbl(iExp, :); 
            newRow = innerjoin(newRow, expMd(:, {'expID', 'expName', 'volumeRate', 'nTrials'}));
            newRow.fl = {allFl'};
            newRow.fwSpeed = {allFwSpeed'};
            newRow.yawSpeed = {allYawSpeed'};
            newRow.odorVols = {allOdorVols'};
            newRow.odorRespFilter = {odorRespFilter'};
            newRow.odorRespMeanSpeed = {meanSpeed'};
            newRow.odorStimDur = dt.subset.offsetTime(1) - ...
                    dt.subset.onsetTime(1);
            newRow.odorRespVector = {odorRespVector'};
            newRow.volTimes = {allVolTimes};
            
            % Append to source data table
            obj.sourceData = [obj.sourceData; newRow];
            
        end%iExp
        
    end% Constructor
    
    % Initialize model data
    function obj = initialize_models(obj, modelParams)
        
        % modelParams = struct with fields: trainTestSplit, criterion, pEnter, pRemove,
        %                                   verbose, useYaw, useDriftCorrection
        % Can be either a scalar structure or have one set of params for each experiment in rm.

        disp('Initializing model data...')
        obj.modelParams = modelParams;
        obj.modelData = [];
        
        for iExp = 1:size(obj.sourceData, 1)
            if numel(modelParams) > 1 
                mp = modelParams(iExp);
            else
                mp = modelParams;
            end
            
            currExpData = obj.sourceData(iExp, :);
            newRow = table(currExpData.expID, 'VariableNames', {'expID'});            
            disp(currExpData.expID{:})
            
            % ------- Create intial table of predictor and response variables --------
            tbl = table(abs(currExpData.fwSpeed{:}'), 'variableNames', {'fwSpeed'});
            if mp.useYaw
                tbl.yawSpeed = currExpData.yawSpeed{:}';
            end
            tbl.odorResp = currExpData.odorRespVector{:}';
            if mp.useDriftCorrection
                tbl.volsFromExpStart = (1:size(tbl, 1))';
            end
            
            % Add fluorescence data in last column of input table
            tbl.fl = currExpData.fl{:}';
            if ~isempty(currExpData.skipVols{:})
                tbl.fl(currExpData.skipVols{:}) = nan;
            end
            
            % Scale all variables so that the max value is either 1 or -1 (and center?)
            for iCol = 1:(size(tbl, 2))
                %         tbl{:, iCol} = tbl{:, iCol} - mean(tbl{:, iCol}, 'omitnan'); % Center mean of all variables at zero
                tbl{:, iCol} = tbl{:, iCol} ./ max(abs(tbl{:, iCol}), [], 'omitnan');
            end
            
            % Add to model data table for future use in plotting model predictions
            newRow.fullDataTbl = {tbl};
            shuffleInds = randperm(size(tbl, 1));
            splitInd = floor(numel(shuffleInds) * mp.trainTestSplit);
            newRow.fitRowInds = {shuffleInds(1:splitInd)};
            newRow.testRowInds = {shuffleInds((splitInd + 1):end)};
            obj.modelData = [obj.modelData; newRow];
            
        end% iExp
        disp('All model data ready for training')
    end% Function
    
    % Plot the trial-averaged odor response filter for a specific experiment
    function ax = plot_odor_filter(obj, expID, axesHandle)
        if nargin < 3
            f = figure(22); clf; hold on;
            f.Color = [1 1 1];
            axesHandle = gca();
        else
            axes(axesHandle);
        end
        if isnumeric(expID)
           expID = obj.sourceData.expID{expID}; 
        end
        currExpData = obj.sourceData(strcmp(obj.sourceData.expID, expID), :);
        odorRespFilter = currExpData.odorRespFilter{:};
        xx = currExpData.volTimes{:}(1:numel(odorRespFilter));
        plot(xx, odorRespFilter, 'linewidth', 3, 'color', 'b');            
        plot_stim_shading([0, currExpData.odorStimDur])
        title([regexprep(expID, '_', '\\_'), ' odor response filter'])
        xlabel('Time (sec)');
        set(gca, 'FontSize', 11)
        ax = axesHandle;
    end
    
    % Plot the trial-averaged moveSpeed response filter for a specific experiment
    function ax = plot_mean_moveSpeed(obj, expID, axesHandle)
        if nargin < 3
            f = figure(23); clf; hold on;
            f.Color = [1 1 1];
            axesHandle = gca();
        else
            axes(axesHandle);
        end
        if isnumeric(expID)
           expID = obj.sourceData.expID{expID}; 
        end
        currExpData = obj.sourceData(strcmp(obj.sourceData.expID, expID), :);
        meanSpeed = currExpData.odorRespMeanSpeed{:};
        frameTimes = (1:1000) ./ 25;
        xx = frameTimes(1:numel(meanSpeed));
        plot(xx, meanSpeed, 'linewidth', 3, 'color', 'k');            
        plot_stim_shading([0, currExpData.odorStimDur])
        title([regexprep(expID, '_', '\\_'), ' mean moveSpeed response'])
        xlabel('Time (sec)');
        set(gca, 'FontSize', 11)
        ax = axesHandle;
    end
    
end%methods


methods(Static)
    
    % Create plot of model coefficients
    function plot_coeffs(mdl, axesHandle)
        
        % Create a new figure and axes if none was provided
        if nargin < 2
            figure(1); clf; hold on
            ax = gca();
        else
            ax = axesHandle; 
        end
        
        % Plot data
        nCoeffs = mdl.NumCoefficients;
        plot(2:nCoeffs, mdl.Coefficients.Estimate(2:end), 'o', 'color', 'k', 'linewidth', 3, ...
            'markersize', 8);
        plot([0 nCoeffs + 1], [0 0], '--', 'color', 'b');
        for iCoeff = 2:nCoeffs
            plot([iCoeff, iCoeff], [0, mdl.Coefficients.Estimate(iCoeff)], '-.', 'color', 'k')
        end
        
        % Format axes
        ax.XLim = [1.5, nCoeffs + 0.5];
        ax.FontSize = 12;
        ax.XTick = 2:numel(mdl.Coefficients.Estimate);
        ax.XTickLabel = regexprep(mdl.CoefficientNames(2:end), '_', '\\_');
        ax.XTickLabelRotation = 30;             
    end
    
end%Static methods

end%class

% ==================================================================================================
% Local functions
% ==================================================================================================
