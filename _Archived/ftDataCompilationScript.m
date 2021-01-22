
% Load behavior experiment

expDates = {'2018_12_03_exp_1', ...
            '2018_12_03_exp_2', ...
            '2018_12_03_exp_3', ...
            '2018_12_04_exp_2', ...
            '2018_12_04_exp_3', ...
            '2018_12_13_exp_1', ...
            '2018_12_13_exp_2', ...
            '2018_12_13_exp_3', ...
            '2018_12_13_exp_4', ...
            '2018_12_16_exp_1', ...
            '2018_12_17_exp_1', ...
            '2018_12_17_exp_2', ...
            '2018_12_17_exp_3', ...
            '2018_12_18_exp_1', ...
            '2018_12_18_exp_2', ...
            '2018_12_18_exp_3', ...
            '2019_01_07_exp_1', ...
            '2019_01_09_exp_1', ...
            '2019_01_09_exp_2', ...
            '2019_01_10_exp_1', ...
            '2019_01_10_exp_2', ...
            '2019_01_11_exp_1', ...
            '2019_01_11_exp_2', ...
    };
sids = zeros(1, numel(expDates));






allExpData = [];
for iExp = 1:numel(expDates)
    
    expDate = expDates{iExp};
    disp(expDate);
    sid = sids(iExp);
    parentDir = fullfile('D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate);
    savePath = fullfile('D:\Dropbox (HMS)\2P Data\Behavior Vids\', expDate, ['sid_', num2str(sid)]);
    
    annotFileName = '_Movies\autoAnnotations.mat';
    
    % ----------------  Load stim metadata -------------------------------------------------------------
    infoStruct = [];
    stimDataFiles = dir(fullfile(parentDir, ['metadata*sid_', num2str(sid), '*.mat']));
    infoStruct.stimOnsetTimes = []; infoStruct.stimDurs = []; infoStruct.trialType = []; infoStruct.outputData = [];
    for iFile = 1:numel(stimDataFiles)
        
        % Load file
        load(fullfile(parentDir, stimDataFiles(iFile).name)); % variable "metaData" with fields 'trialDuration', 'nTrials', 'stimTypes', 'sid', 'taskFile', 'outputData'
        
        % Add block data
        infoStruct.blockData(iFile).nTrials = metaData.nTrials;
        infoStruct.blockData(iFile).stimTypes = metaData.stimTypes;
        infoStruct.blockData(iFile).taskFile = metaData.taskFile;
%         infoStruct.blockData(iFile).outputData = metaData.outputData;
        
        % Add info to overall session data
        infoStruct.trialDuration = metaData.trialDuration;
        infoStruct.expDate = expDate;
        
        % Add downsampled output data broken down into trials
        currOutput = metaData.outputData';
        sampPerTrial = size(currOutput, 2) / metaData.nTrials;
        rsOutput = reshape(currOutput, size(currOutput, 1), sampPerTrial, metaData.nTrials);   % --> [channel, sample, trial]
        disp(size(rsOutput))
        infoStruct.outputData = cat(3, infoStruct.outputData, permute(rsOutput(:,1:100:end,:), [2 1 3]));    % --> [sample, channel, trial]
        
        % Add trial type info
        for iTrial = 1:metaData.nTrials
            currStim = metaData.stimTypes{iTrial};
            infoStruct.stimOnsetTimes(end + 1) = str2double(regexp(currStim, '(?<=Onset_).*(?=_Dur)', 'match'));
            infoStruct.stimDurs(end + 1) = str2double(regexp(currStim, '(?<=Dur_).*', 'match'));
            infoStruct.trialType{end + 1} = regexp(currStim, '.*(?=_Onset)', 'match', 'once');
            if strcmp(infoStruct.trialType{end}, 'NoOdor')
                infoStruct.trialType{end} = 'NoStim'; % For backwards compatibility
            end
        end%iTrial
    end%iFile
    infoStruct.nTrials = size(infoStruct.outputData, 3);
    infoStruct.stimTypes = sort(unique(infoStruct.trialType));
    infoStruct.stimSepTrials = [];
    for iStim = 1:length(infoStruct.stimTypes)
        infoStruct.stimSepTrials.(infoStruct.stimTypes{iStim}) = logical(cellfun(@(x) ...
            strcmp(x, infoStruct.stimTypes{iStim}), infoStruct.trialType));
    end
    
    % ----------------  Load autoAnnotation data -------------------------------------------------------
    annotData = load(fullfile(parentDir, annotFileName)); % variables 'trialAnnotations', 'annotParams', 'ftData', 'flowArr', 'goodTrials', 'behaviorLabels', 'frameInfo'
    annotData.nFrames = annotData.frameInfo.nFrames;
    annotData.frameTimes = annotData.frameInfo.frameTimes;
    annotData.FRAME_RATE = annotData.frameInfo.FRAME_RATE;
    infoStruct = setstructfields(infoStruct, annotData);
    
    
    % ----------------  Load FicTrac data --------------------------------------------------------------
    ftData = load_fictrac_data(infoStruct, 'sid', sid, 'ParentDir', fullfile(parentDir, '\_Movies\FicTracData'));
    infoStruct.ftData = ftData;
    infoStruct.goodTrials(logical(ftData.resets)) = 0;
    
    % ---------------- Create workspace vars -----------------------------------------------------------
    infoStruct = orderfields(infoStruct);
    
    clear annotData currOutput rsOutput metaData ftData
    allExpData{iExp} = infoStruct;
end


%%

save(fullfile('D:\Dropbox (HMS)\2P Data\Behavior Vids', 'allBehaviorData.mat'), 'expDates', 'allExpData', '-v7.3')

%%
load(fullfile('D:\Dropbox (HMS)\2P Data\Behavior Vids', 'allBehaviorData.mat'), 'expDates', 'allExpData')
%%

FRAME_RATE = 25;

allMeanSpeeds = []; medSpeeds = []; totalAvgSpeeds = [];
for iExp = 1:numel(allExpData)
   
    currExpData = allExpData{iExp};
    currFtData = currExpData.ftData.moveSpeed * FRAME_RATE * 4.5; % Convert to mm/sec
    
    meanSpeed = movmean(squeeze(mean(currFtData, 1)), 3);
    
    allMeanSpeeds{iExp} = meanSpeed;
    medSpeeds(iExp) = median(meanSpeed);
    totalAvgSpeeds(iExp) = mean(meanSpeed);
end




% figure(1); clf; hold on
% for iExp = 1:numel(allMeanSpeeds)
%    plot(allMeanSpeeds{iExp});
%    
%     
% end
% legend(cellfun(@regexprep, expDates, repmat({'_'}, 1, numel(expDates)), repmat({'\\_'}, 1, numel(expDates)), 'UniformOutput', false));








