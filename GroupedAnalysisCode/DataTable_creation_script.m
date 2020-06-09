%% Load a base aligned data object

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\';
alignEventDateStr = '20200608';



% % ballstop, flailing, grooming, isolatedmovement, locomotion, odor, optostim, panelsflash, soundstim
% alignEventName = 'panelsflash';

eventNames = {'ballstop', 'flailing', 'grooming', 'isolatedmovement', 'locomotion', 'odor', ...
        'optostim', 'panelsflash', 'soundstim'};

for iEvent = 1:numel(eventNames)

    alignEventName = eventNames{iEvent};
    
    load(fullfile(parentDir, 'Saved_AlignEvent_objects', [alignEventDateStr, '_AlignEventObj_', ...
        alignEventName, '.mat']));
    
    % Save data tables with several different common alignment windows
    preStimWinSizes = [2];
    postStimWinSizes = [1 2 3 4];
    
    for iPre = 1:numel(preStimWinSizes)
        for iPost = 1:numel(postStimWinSizes)
            analysisWin = [preStimWinSizes(iPre), postStimWinSizes(iPost)];
            dt = alignObj.output_analysis_DataTable(analysisWin);
            saveDir = fullfile(parentDir, 'Saved_DataTable_objects');
            fileName = [alignEventName, '_pre_', num2str(analysisWin(1)), '_post_', ...
                num2str(analysisWin(2)), '.mat'];
            save(fullfile(saveDir, fileName), 'dt');
        end
    end
 
end