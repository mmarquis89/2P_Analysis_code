


%%
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';

expList = load_expList();

%%

iExp = 137;

currExpID = expList.expID{iExp};
disp(currExpID)
roiFile = fullfile(parentDir, [currExpID, '_roiData.mat']);
trialMdFile = fullfile(parentDir, [currExpID, '_trialMetadata.csv']);
if exist(trialMdFile, 'file')
    trialMd = readtable(trialMdFile, 'delimiter', ',');
    trialMd.pmtShutoffVols = num2cell(nan(size(trialMd, 1), 1));
    if exist(roiFile, 'file')
        
        % Load ROI data for current experiment
        load(roiFile, 'roiData');
        
        % Get names of all ROIs
        uniqueNames = unique(roiData.roiName);
        
        % Average data from all ROIs together, concatenating all trials
        for iName = 1:numel(uniqueNames)
            currRoiData = roiData(strcmp(roiData.roiName, uniqueNames{iName}), :);
            if iName == 1
                roiDataVect = cell2mat(currRoiData.rawFl');
                
                % Save trial bounds for later
                trialStartVols = 1 + cumsum(cellfun(@numel, currRoiData.rawFl));
                trialStartVols = [1, trialStartVols(1:end-1)'];
                
            else
                if numel(cell2mat(currRoiData.rawFl)) == numel(roiDataVect)
                    roiDataVect = roiDataVect + cell2mat(currRoiData.rawFl');
                end
            end
        end
        
        % Plot volume-to-volume differences in averaged ROI data
        figure(1); clf; hold on
        plot(diff(roiDataVect ./ numel(uniqueNames)));
        xlim([0 numel(roiDataVect)])
        
        % Request onset volumes of PMT shutoffs from user
        offVols = str2num(input('Enter all PMT shutoff volumes, separated by spaces:\n', 's'));
        
        if ~isempty(offVols)
            
            % Request offset volumes
            onVols = str2num(input('Enter all PMT turn-on volumes, separated by spaces:\n', 's'));
            
            % Create logical vector of PMT shutoffs
            expShutoffVols = zeros(size(roiDataVect));
            for iVol = 1:numel(offVols)
                expShutoffVols(offVols(iVol):onVols(iVol)) = 1;
            end
            
            % Split into trials and save to trialMetadata table
            trialBounds = [trialStartVols, numel(roiDataVect)];
            for iTrial = 1:numel(trialStartVols)
                trialMd.pmtShutoffVols{iTrial} = ...
                        expShutoffVols(trialBounds(iTrial):trialBounds(iTrial + 1));
            end
            
        end
    end%if
    
    % Re-save trialMetadata table
    writetable(trialMd, trialMdFile);
    
    disp(['PMT shutoff processing completed for ', currExpID])
end%if
