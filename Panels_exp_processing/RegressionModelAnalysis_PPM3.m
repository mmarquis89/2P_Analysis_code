classdef RegressionModelAnalysis_PPM3
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
    function obj = RegressionModelAnalysis_PPM3(expInfoTbl, sourceDataParams)
        
        % expInfoTable = table with colums: [expID][skipTrials][skipVols]
        
        % sourceDataParams = struct with fields: maxSpeed, smWinVols, smWinFrames, smReps, 
        %                                        ftLagVols, speedType
        
        obj.sourceDataParams = sourceDataParams;
        
        % Unpack sourceDataParams
        maxSpeed = sourceDataParams.maxSpeed;
        smWinVols = sourceDataParams.smWinVols;
        smWinFrames = sourceDataParams.smWinFrames;
        smReps = sourceDataParams.smReps;
        ftLagVols = sourceDataParams.ftLagVols;
        roiName = sourceDataParams.roiName;
        speedType = sourceDataParams.speedType;
        obj.modelData = [];
        
        % -------  Load source data for all experiments -------
        dataDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\EB-DAN_GroupedAnalysisData';
        FRAME_RATE = 25;
        
        % Pre-process data for each experiment
        obj.sourceData = [];
        for iExp = 1:numel(expInfoTbl.expID)
            
            currExpID = expInfoTbl.expID{iExp};
            
            % -------  Load source data -------
            [expMd, trialMd] = load_metadata({currExpID}, dataDir);
            roiData = load_roi_data({currExpID}, dataDir);
            roiNames = regexp(roiData.roiName, roiName, 'match', 'once'); % Just one hemi
            roiNames(cellfun(@isempty, roiNames)) = [];   
            roiData = roiData(strcmp(roiData.roiName, roiNames{1}), :);
            ft = load_ft_data({currExpID}, dataDir);
            
            % Process data for each trial (assuming no time between trials)
            trialNums = trialMd.trialNum;
            if ~isempty(expInfoTbl.skipTrials{iExp})
                trialNums(ismember(trialNums, expInfoTbl.skipTrials{iExp})) = [];
            end
            allSpeed = []; allYawSpeed = []; allFl = []; allVolTimes = [];
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
                if strcmp(speedType, 'moveSpeed')
                    currSpeed = repeat_smooth(currFt.moveSpeed{:}, smReps, ...
                            'smWin', smWinFrames);
                else
                    currSpeed = repeat_smooth(currFt.fwSpeed{:}, smReps, ...
                            'smWin', smWinFrames);
                end
                currSpeed(currSpeed > maxSpeed) = maxSpeed;
                currYawSpeed = abs(repeat_smooth(rad2deg(currFt.yawSpeed{:}), ...
                        smReps, 'smWin', smWinFrames));
                currSpeedVols = []; currYawSpeedVols = [];
                for iVol = 1:numel(volTimes)
                    dsFrame = argmin(abs(currFt.frameTimes{:} - volTimes(iVol)));
                    currSpeedVols(iVol) = currSpeed(dsFrame);
                    currYawSpeedVols(iVol) = currYawSpeed(dsFrame);
                end
                
                % Apply a lag to the FicTrac data
                currSpeedVols = currSpeedVols([ones(1, ftLagVols), 1:end-ftLagVols]);
                currYawSpeedVols = currYawSpeedVols([ones(1, ftLagVols), 1:end-ftLagVols]);
                
                % Append variables to whole-experiment vectors
                allFl = [allFl; currFl];
                allSpeed = [allSpeed; currSpeedVols'];
                allYawSpeed = [allYawSpeed; currYawSpeedVols'];
                if iTrial == 1
                    allVolTimes = volTimes;
                else
                    allVolTimes = [allVolTimes, volTimes + allVolTimes(end)];
                end
            end%iTrial
            
            % Add current experiment to sourceData table
            newRow = expInfoTbl(iExp, :); 
            newRow = innerjoin(newRow, expMd(:, {'expID', 'expName', 'volumeRate', 'nTrials'}));
            newRow.fl = {allFl'};
            newRow.(speedType) = {allSpeed'};
            newRow.yawSpeed = {allYawSpeed'};
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
