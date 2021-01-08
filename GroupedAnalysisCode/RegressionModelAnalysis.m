classdef RegressionModelAnalysis
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
%       ax = plot_odor_filter(expID, axesHandle)       
%       plot_mean_moveSpeed(expID, axesHandle)
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
    function obj = RegressionModelAnalysis(expInfoTbl, sourceDataParams)
        
        % expInfoTable = table with colums: [expID][skipTrials][skipVols]
        
        % sourceDataParams = struct with fields: roiName 
        %                                        maxSpeed 
        %                                        smWinVols
        %                                        smWinFrames
        %                                        smReps 
        %                                        ftLagVols 
        %                                        speedType 
        %                                        speedType
        %                                        odorRespFilterDur
        %                                        parentDir
        %                                        dataDir
        %                                        eventDataParentDir
        %                                        alignEventDateStr
        %                                        convertFtUnits
        %                                        loadOneNoteData
        %                                        alignObjFilterDefs
        
        p = sourceDataParams;
        
        obj.modelData = [];
        
        % -------  Load source data for all experiments -------
        if ~isfield(p, 'parentDir') || isempty(p.parentDir)
            p.parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData';
        end
        if ~isfield(p, 'dataDir') || isempty(p.dataDir)
            p.dataDir = fullfile(p.parentDir, 'all_experiments');
        end
                
        % Load AlignEvent object
        disp('Loading new AlignEvent object...');
        if ~isfield(p, 'eventDataParentDir') || isempty(p.eventDataParentDir)
            p.eventDataParentDir = p.parentDir;
        end
        if ~isfield(p, 'alignEventDateStr') || isempty(p.alignEventDateStr)
            p.alignEventDateStr = '20200609';
        end
        load(fullfile(p.eventDataParentDir, 'Saved_AlignEvent_objects', [p.alignEventDateStr, ...
                '_AlignEventObj_odor.mat']), 'alignObj');
        disp('Loading complete')
        
        % Load event DataTable
        disp('Loading event DataTable...')
        fileName = ['odor_pre_', num2str(sourceDataParams.odorRespFilterDur(1)), '_post_', ...
                num2str(sourceDataParams.odorRespFilterDur(2)), '.mat'];
        load(fullfile(p.eventDataParentDir, 'Saved_DataTable_objects', fileName), 'dt')
        disp('Loading complete')
        
        % Pre-process data for each experiment
        obj.sourceData = [];
        for iExp = 1:numel(expInfoTbl.expID)
            
            currExpID = expInfoTbl.expID{iExp};
            
            % -------  Load source data -------
            [expMd, trialMd] = load_metadata({currExpID}, p.dataDir);
            roiData = load_roi_data({currExpID}, p.dataDir);
            roiNames = regexp(roiData.roiName, [p.roiName, '-.'], 'match', 'once');
            roiNames(cellfun(@isempty, roiNames)) = [];   
            roiData = roiData(strcmp(roiData.roiName, roiNames{1}), :);
            ft = load_ft_data({currExpID}, p.dataDir);
            
            % Load odor event data
            eventData = load_event_data({currExpID}, p.dataDir);
            odorEventData = eventData.odor;
            
            % Process data for each trial (assuming no time between trials)
            trialNums = trialMd.trialNum;
            if ~isempty(expInfoTbl.skipTrials{iExp})
                trialNums(ismember(trialNums, expInfoTbl.skipTrials{iExp})) = [];
            end
            allSpeed = []; allYawSpeed = []; allOdorVols = []; allFl = []; allVolTimes = [];
            for iTrial = 1:numel(trialNums)
                currTrialNum = trialNums(iTrial);
                md = innerjoin(expMd, trialMd(currTrialNum, :));
                volTimes = calc_volTimes(md.nVolumes, md.volumeRate, md.trialDuration, ...
                        md.originalTrialCount);
                
                % Smooth imaging rawFl data for current trial
                currRoiData = roiData(currTrialNum, :);
                currFl = smoothdata(currRoiData.rawFl{:}, 'gaussian', p.smWinVols);
                if ~isnan(md.pmtShutoffVols{:})
                    currFl(md.pmtShutoffVols{:}) = nan;
                end
                
                % Smooth FicTrac data and convert units if necessary (for older experiments)
                currFt = ft(ft.trialNum == currTrialNum, :);
                if strcmp(p.speedType, 'moveSpeed')
                    currSpeed = repeat_smooth(currFt.moveSpeed{:}, p.smReps, ...
                            'smWin', p.smWinFrames);
                else
                    currSpeed = repeat_smooth(currFt.fwSpeed{:}, p.smReps, ...
                            'smWin', p.smWinFrames);
                end
                if p.convertFtUnits
                    currSpeed = currSpeed .* 4.5 .* 25; % Convert units if necessary
                end
                currSpeed(currSpeed > p.maxSpeed) = p.maxSpeed;
                currYawSpeed = abs(repeat_smooth(currFt.yawSpeed{:}, p.smReps, 'smWin', ...
                        p.smWinFrames));
                if p.convertFtUnits
                    currYawSpeed = rad2deg(currYawSpeed) * 25; % Convert units if necessary
                end
                
                % Downsample to match volume rate
                currSpeedVols = []; currYawSpeedVols = [];
                for iVol = 1:numel(volTimes)
                    dsFrame = argmin(abs(currFt.frameTimes{:} - volTimes(iVol)));
                    currSpeedVols(iVol) = currSpeed(dsFrame);
                    currYawSpeedVols(iVol) = currYawSpeed(dsFrame);
                end
                
                % Apply a lag to the FicTrac data
                currSpeedVols = currSpeedVols([ones(1, p.ftLagVols), 1:end-p.ftLagVols]);
                currYawSpeedVols = currYawSpeedVols([ones(1, p.ftLagVols), 1:end-p.ftLagVols]);
                
                % Get odor command vector for current trial
                currTrialOdorVols = odorEventData.create_logical_array(md.nVolumes, volTimes, ...
                        table({currExpID}, currTrialNum, 'variableNames', {'expID', 'trialNum'}));
                
                % Append variables to whole-experiment vectors
                allFl = [allFl; currFl];
                allSpeed = [allSpeed; currSpeedVols'];
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
            if ~isfield(p, 'loadOneNoteData') || isempty(p.loadOneNoteData)
                p.loadOneNoteData = 0;
            end
            if ~isfield(p, 'alignObjFilterDefs') || isempty(p.alignObjFilterDefs)
                filterDefs = alignObj.create_filterDefs('loadOneNoteData', p.loadOneNoteData);
                filterDefs.expID = currExpID;
                filterDefs.roiName = roiNames{1};
                filterDefs.odor = -1;
                filterDefs.locomotion = [];
                filterDefs.trialNum = trialNums;
                p.alignObjFilterDefs = filterDefs;
            else
                p.alignObjFilterDefs.expID = currExpID;
                p.alignObjFilterDefs.trialNum = trialNums;
                p.alignObjFilterDefs.roiName = roiNames{1};
            end
            dt = dt.initialize_filters(p.alignObjFilterDefs);

            meanTrace = smoothdata(mean(cell2mat(dt.subset.eventFl'), 2, 'omitnan'), ...
                    'gaussian', p.smWinVols);
            odorRespFilter = meanTrace(dt.subset.volTimes{1} > 0); % Trim pre-onset volumes
            odorRespFilter = odorRespFilter - odorRespFilter(1);   % Offset to zero first volume
            
            % Also calculate a trial-averaged version of the move speed while we're at it
            if strcmp(p.alignObjFilterDefs.roiName, 'TypeF-R')
                filterDefs = p.alignObjFilterDefs;
                filterDefs.expID = currExpID;
                filterDefs.roiName = 'TypeF';
                filterDefs.locomotion = [];
                dt2 = dt.initialize_filters(filterDefs);
            else
                dt2 = dt;
            end
            moveSpeedArr = cell2mat(dt2.subset.moveSpeed');
            if p.convertFtUnits
                moveSpeedArr = moveSpeedArr .* 4.5 .* 25;
            end
            meanSpeed = repeat_smooth(mean(moveSpeedArr, 2, 'omitnan'), p.smReps, 'smWin', ...
                    p.smWinFrames);
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
            newRow.(p.speedType) = {allSpeed'};
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
        
        % Save source data params to object
        obj.sourceDataParams = p;
        
    end% Constructor
    
    % Initialize model data
    function obj = initialize_models(obj, modelParams)
        
        % modelParams = struct with fields: trainTestSplit, 
        %                                   kFold
        %                                   criterion
        %                                   upper
        %                                   pEnter 
        %                                   pRemove
        %                                   verbose
        %                                   useYaw
        %                                   useDriftCorrection 
        %                                   odorIntegrationWin
        %                                   speedPadDist
        %                                   speedIntegrationWin
        %                                   standardizeInputs
        %                                   normalizeInputs
        % Can be either a scalar structure or have one set of params for each experiment in rm.

        disp('Initializing model data...')
        obj.modelParams = modelParams;
        obj.modelData = [];
        speedType = obj.sourceDataParams.speedType;
        
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
            tbl = table(abs(currExpData.(speedType){:}'), 'variableNames', {speedType});
            if mp.useYaw
                tbl.yawSpeed = currExpData.yawSpeed{:}';
            end
            tbl.odorResp = currExpData.odorRespVector{:}';
            if mp.useDriftCorrection
                tbl.volsFromExpStart = (1:size(tbl, 1))';
            end
            
            % Calculate integrated odor values if necessary
            if ~isempty(mp.odorIntegrationWin)
                odorIntegrationWinVols = round(mp.odorIntegrationWin * currExpData.volumeRate);
                integratedOdor = [];
                for iWin = 1:numel(odorIntegrationWinVols)
                    currIntegrationWin = odorIntegrationWinVols(iWin);
                    for iVol = 1:numel(currExpData.odorVols{:})
                        if iVol <= currIntegrationWin
                            integratedOdor(iWin, iVol) = sum(currExpData.odorVols{:}(1:iVol));
                        else
                            integratedOdor(iWin, iVol) = sum(currExpData.odorVols{:}(iVol - ...
                                currIntegrationWin:iVol));
                        end
                    end
                end
                for iWin = 1:size(integratedOdor, 1)
                    varName = ['odorHistory_', num2str(mp.odorIntegrationWin(iWin))];
                    tbl.(varName) = integratedOdor(iWin, :)';
                end
            end
            
            % Calculate mean speed history values if necessary
            padDist = mp.speedPadDist;
            padDistVols = padDist * currExpData.volumeRate;
            if ~isempty(mp.speedIntegrationWin)
                speedIntegrationWinVols = round(mp.speedIntegrationWin * currExpData.volumeRate);
                meanSpeed = [];
                for iWin = 1:numel(speedIntegrationWinVols)
                    currIntegrationWin = speedIntegrationWinVols(iWin);
                    for iVol = 1:numel(currExpData.(speedType){:})
                        if iVol <= padDistVols
                            meanSpeed(iWin, iVol) = sum(currExpData.(speedType){:}(1:iVol));
                        elseif iVol <= currIntegrationWin + padDistVols
                            meanSpeed(iWin, iVol) = sum(currExpData.(speedType){:}(1:(iVol - ...
                                    padDistVols)));
                        else
                            meanSpeed(iWin, iVol) = sum(currExpData.(speedType){:}(iVol - ...
                                currIntegrationWin:(iVol - padDistVols)));
                        end
                    end
                end
                for iWin = 1:size(meanSpeed, 1)
                    varName = ['speedHistory_', num2str(mp.speedIntegrationWin(iWin))];
                    tbl.(varName) = meanSpeed(iWin, :)';
                end
            end
            
            % Add fluorescence data in last column of input table
            tbl.fl = currExpData.fl{:}';
            if ~isempty(currExpData.skipVols{:})
                tbl.fl(currExpData.skipVols{:}) = nan;
            end
            
            % Standardize or normalize the input variables as needed
            for iCol = 1:(size(tbl, 2))
                
                if mp.standardizeInputs
                    if mp.normalizeInputs
                        error('"standardizeInputs" and "normalizeInputs" params cannot both be 1"');
                    end
                    
                    % Standardize and center variables by Z-scoring
                    tbl{:, iCol} = normalize(tbl{:, iCol}); 
                    
                elseif mp.normalizeInputs
                    
                    % Normalize inputs by re-scaling from 0-1
                    tbl{:, iCol} = tbl{:, iCol} ./ max(abs(tbl{:, iCol}), [], 'omitnan');  
                end
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
    
    % Find optimal odor integration window for each experiment
    function obj = optimize_odor_integration_windows(obj)
        disp('Finding best odor integration window sizes...')
        
        bestWinSizes = [];
        cvModels = {};
        cvAdjR2Vals = {};
        for iExp = 1:size(obj.sourceData, 1)
            if numel(obj.modelParams) > 1
                mp = obj.modelParams(iExp);
            else
                mp = obj.modelParams;
            end
            disp(obj.sourceData.expID{iExp})
            
            currExpData = obj.modelData(iExp, :);
            
            % Get training data for current experiment
            tblTrain = currExpData.fullDataTbl{:}(currExpData.fitRowInds{:}, :);
            
            % Drop rows without valid fluorescence measurements
            tblTrain = tblTrain(~isnan(tblTrain.fl), :);
            
            % Drop yawSpeed variable
%             tblTrain = tblTrain(:, ~strcmp(tblTrain.Properties.VariableNames, 'yawSpeed'));
            
            % Divide training data into cross-validation folds
            chunkSize = floor(size(tblTrain, 1) / mp.kFold);
            startInds = 1:chunkSize:size(tblTrain, 1);
            endInds = [startInds(2:end) - 1, size(tblTrain, 1)];
            
            % Identify positions of odor history variables
            kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
                mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
            emptyArgs = cellfun(@isempty, kvArgs);
            kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
            varNames = tblTrain.Properties.VariableNames;
            odorHistVars = ~cellfun(@isempty,regexp(varNames, 'odorHistory'));
            odorHistVarInds = find(odorHistVars);
            
            % Select odor history window size using cross-validation
            nWindows = numel(mp.odorIntegrationWin);
            loopTbls = [];
            for iWin = 1:nWindows
                loopTbls{iWin} = tblTrain(:, [1:(odorHistVarInds(1) - 1), odorHistVarInds(iWin), ...
                    size(tblTrain, 2)]);
            end
            kFold = mp.kFold;
            predAdjR2 = zeros(nWindows, kFold);
            allCvMdls = cell(nWindows, kFold);
            parfor iWin = 1:nWindows

                currTbl = loopTbls{iWin};
                
                % Fit model with cross-validation
                for iFold = 1:kFold
                    if iWin == 1 && ~mod(iFold, 20)
                        disp(['Running cross-validation fold ', num2str(iFold), ' of ', ...
                                num2str(kFold)]);
                    end
                    
                    % Split data for current fold
                    currFoldTrainTbl = currTbl;
                    currFoldTrainTbl(startInds(iFold):endInds(iFold), :) = [];
                    currFoldTestTbl = currTbl(startInds(iFold):endInds(iFold), :);
                    
                    % ------------- Generate model ---------------
                    mdl = stepwiselm(currFoldTrainTbl, kvArgs{:});
                    % mdl = fitlm(tblFit);
                    allCvMdls{iWin, iFold} = compact(mdl);
                    
                    % Evaluate fit to test data
                    [yPred, ~] = predict(mdl, currFoldTestTbl(:, 1:end-1));
                    [~, R2] = r_squared(currFoldTestTbl.fl, yPred, mdl.NumCoefficients);
                    predAdjR2(iWin, iFold) = R2;
                    
                end%iFold
            end%iWin
            
            % Store result variables
            cvModels{iExp, 1} = allCvMdls;
            cvAdjR2Vals{iExp, 1} = predAdjR2;
            bestWinSizes(iExp, 1) = mp.odorIntegrationWin(argmax(mean(predAdjR2, 2), 1));
            
        end%iExp
        
        % Save results for all experiments to model data table
        obj.modelData.cvModels = cvModels;
        obj.modelData.cvAdjR2Vals = cvAdjR2Vals;
        obj.modelData.bestWinSizes = bestWinSizes;        
    
        disp('Optimal odor integration windows selected');
        
    end% function
 
    % Find optimal mean speed window for each experiment
    function obj = optimize_mean_speed_windows(obj)
        disp('Finding best mean speed window sizes...')
        
        bestWinSizes = [];
        cvModels = {};
        cvAdjR2Vals = {};
        for iExp = 1:size(obj.sourceData, 1)
            if numel(obj.modelParams) > 1
                mp = obj.modelParams(iExp);
            else
                mp = obj.modelParams;
            end
            disp(obj.sourceData.expID{iExp})
            
            currExpData = obj.modelData(iExp, :);
            
            % Get training data for current experiment
            tblTrain = currExpData.fullDataTbl{:}(currExpData.fitRowInds{:}, :);
            
            % Drop rows without valid fluorescence measurements
            tblTrain = tblTrain(~isnan(tblTrain.fl), :);
            
            % Drop yawSpeed variable
%             tblTrain = tblTrain(:, ~strcmp(tblTrain.Properties.VariableNames, 'yawSpeed'));
            
            % Divide training data into cross-validation folds
            chunkSize = floor(size(tblTrain, 1) / mp.kFold);
            startInds = 1:chunkSize:size(tblTrain, 1);
            endInds = [startInds(2:end) - 1, size(tblTrain, 1)];
            
            % Identify positions of mean speed history variables
            kvArgs = {'criterion', mp.criterion, 'pEnter', mp.pEnter, 'pRemove', ...
                mp.pRemove, 'verbose', mp.verbose, 'upper', mp.upper};
            emptyArgs = cellfun(@isempty, kvArgs);
            kvArgs(logical(emptyArgs + [emptyArgs(2:end), 0])) = [];
            varNames = tblTrain.Properties.VariableNames;
            speedHistVars = ~cellfun(@isempty,regexp(varNames, 'speedHistory'));
            speedHistVarInds = find(speedHistVars);
            
            % Select speed history window size using cross-validation
            nWindows = numel(mp.speedIntegrationWin);
            loopTbls = [];
            for iWin = 1:nWindows
                loopTbls{iWin} = tblTrain(:, [1:(speedHistVarInds(1) - 1), speedHistVarInds(iWin), ...
                    size(tblTrain, 2)]);
            end
            kFold = mp.kFold;
            predAdjR2 = zeros(nWindows, kFold);
            allCvMdls = cell(nWindows, kFold);
            parfor iWin = 1:nWindows

                currTbl = loopTbls{iWin};
                
                % Fit model with cross-validation
                for iFold = 1:kFold
                    if iWin == 1 && ~mod(iFold, 20)
                        disp(['Running cross-validation fold ', num2str(iFold), ' of ', ...
                                num2str(kFold)]);
                    end
                    
                    % Split data for current fold
                    currFoldTrainTbl = currTbl;
                    currFoldTrainTbl(startInds(iFold):endInds(iFold), :) = [];
                    currFoldTestTbl = currTbl(startInds(iFold):endInds(iFold), :);
                    
                    % ------------- Generate model ---------------
                    mdl = stepwiselm(currFoldTrainTbl, kvArgs{:});
                    % mdl = fitlm(tblFit);
                    allCvMdls{iWin, iFold} = compact(mdl);
                    
                    % Evaluate fit to test data
                    [yPred, ~] = predict(mdl, currFoldTestTbl(:, 1:end-1));
                    [~, R2] = r_squared(currFoldTestTbl.fl, yPred, mdl.NumCoefficients);
                    predAdjR2(iWin, iFold) = R2;
                    
                end%iFold
            end%iWin
            
            % Store result variables
            cvModels{iExp, 1} = allCvMdls;
            cvAdjR2Vals{iExp, 1} = predAdjR2;
            bestWinSizes(iExp, 1) = mp.speedIntegrationWin(argmax(mean(predAdjR2, 2), 1));
            
        end%iExp
        
        % Save results for all experiments to model data table
        obj.modelData.cvModels = cvModels;
        obj.modelData.cvAdjR2Vals = cvAdjR2Vals;
        obj.modelData.bestWinSizes = bestWinSizes;        
    
        disp('Optimal speed integration windows selected');
        
    end% function
    
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

