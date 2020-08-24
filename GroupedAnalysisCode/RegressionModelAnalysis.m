classdef RegressionModelAnalysis
% ==================================================================================================
%
% Properties:
%
%
% Methods:
%       initialize_models()
%       optimize_integration_windows()
%       
%       plot_odor_filter();       
%
% ==================================================================================================

properties
    sourceData
    sourceDataParams
    modelData
    modelParams
end


properties (Dependent)
    
end 


methods
    
    % Constructor
    function obj = RegressionModelAnalysis(expInfoTbl, sourceDataParams)
        
        % expInfoTable = table with colums: [expID][skipTrials][skipVols]
        
        % sourceDataParams = struct with fields: maxSpeed, smWinVols, smWinFrames, smReps, 
        %                                        ftLagVols, odorRespFilterDur
        
        obj.sourceDataParams = sourceDataParams;
        
        % Unpack sourceDataParams
        maxSpeed = sourceDataParams.maxSpeed;
        smWinVols = sourceDataParams.smWinVols;
        smWinFrames = sourceDataParams.smWinFrames;
        smReps = sourceDataParams.smReps;
        ftLagVols = sourceDataParams.ftLagVols;
        
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
            roiNames = regexp(roiData.roiName, ['TypeD', '-.'], 'match', 'once');
            roiNames(cellfun(@isempty, roiNames)) = [];   
            roiData = roiData(strcmp(roiData.roiName, roiNames{1}), :);
            ft = load_ft_data({currExpID}, dataDir);
            
            % Load odor event data
            eventData = load_event_data({currExpID}, dataDir);
            odorEventData = eventData.odor;
            
            % Process data for each trial (assuming no time between trials)
            trialNums = trialMd.trialNum;
            if ~isempty(expInfoTbl.skipTrials{iExp})
                trialNums(trialNums == skipTrials) = [];
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
                currFwSpeed = repeat_smooth(currFt.fwSpeed{:} .* 4.5 .* FRAME_RATE, smReps, 'smWin', ...
                    smWinFrames);
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
            filterDefs.expID = currExpID;
            filterDefs.roiName = roiNames{1};
            filterDefs.odor = -1;
            filterDefs.locomotion = 0;
            filterDefs.flailing = 0;
            filterDefs.grooming = 0;
            dt = dt.initialize_filters(filterDefs);
            if size(dt.subset, 1) == 0
                % Relax the locomotion exclusion requirement if there isn't any data for this experiment
                filterDefs.locomotion = [];
                dt = dt.initialize_filters(filterDefs);
                disp('No movement free odor stim windows - using full exp trial-averaged response')
            else
                disp(['Using ', num2str(size(dt.subset, 1)), ...
                        ' trials to estimate odor filter'])
            end
            meanTrace = smoothdata(mean(cell2mat(dt.subset.eventFl'), 2, 'omitnan'), ...
                    'gaussian', smWinVols);
            odorRespFilter = meanTrace(dt.subset.volTimes{1} > 0); % Trim pre-onset volumes
            odorRespFilter = odorRespFilter - odorRespFilter(1);          % Offset to zero first volume
            
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
        
        % modelParams = struct with fields: trainTestSplit, kFold, criterion, pEnter, pRemove,
        %                                   verbose, odorIntegrationWin
        % Can be either a scalar structure or have one set of params for each experiment in rm.

        disp('Initializing model data...')
        obj.modelParams = modelParams;
        
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
            tbl.yawSpeed = currExpData.yawSpeed{:}';
            tbl.odorResp = currExpData.odorRespVector{:}';
            tbl.volsFromExpStart = (1:size(tbl, 1))';
            
            % Calculate integrated odor values
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
            
            % Add fluorescence data in last column of input table
            tbl.fl = currExpData.fl{:}';
            if ~isempty(currExpData.skipVols{:})
                tbl.fl(currExpData.skipVols{:}) = nan;
            end
            
            % Scale all variables so that the max value is 1 (and center?)
            for iCol = 1:(size(tbl, 2))
                %         tbl{:, iCol} = tbl{:, iCol} - mean(tbl{:, iCol}, 'omitnan'); % Center mean of all variables at zero
                tbl{:, iCol} = tbl{:, iCol} ./ max(tbl{:, iCol}, [], 'omitnan');
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
    function obj = optimize_integration_windows(obj)
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
        plot(xx, odorRespFilter, 'linewidth', 3, 'color', 'k');            
        plot_stim_shading([0, currExpData.odorStimDur])
        title([regexprep(expID, '_', '\\_'), ' odor response filter'])
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
