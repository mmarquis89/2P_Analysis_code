function vid_dir_cleanup(vidDataDir, sid)


try
    
    % Update to reflect any changes made to the imaging data during scanimage file cleanup
    volCountFile = ['sid_', num2str(sid), '_volume_counts.mat'];
    if exist(fullfile(vidDataDir, volCountFile), 'file')
        load(fullfile(vidDataDir, volCountFile)); % 'volCounts', 'trialFiles', 'bidList'
        
        
        for iBlock = 1:numel(bidList)
            
            bid = bidList(iBlock);
            disp(['Processing block #', num2str(bid)])
            write_to_log(['Processing block #', num2str(bid)], mfilename);
            trialVidDirs = dir(fullfile(vidDataDir, ['*sid_', num2str(sid), '*bid_', num2str(bid), ...
                '*tid*']));
            trialVidDirs(~[trialVidDirs.isdir]) = [];
            if ~isdir(fullfile(vidDataDir, 'BadTrialsBackup'))
                mkdir(fullfile(vidDataDir, 'BadTrialsBackup'))
            end
            if isempty(dir(fullfile(vidDataDir, 'BadTrialsBackup', ['*_bid_', num2str(bid), '.mat'])))
                
                changeLogFile = ['sid_', num2str(sid), '_bid_', num2str(bid), ...
                    '_single_vol_trials.mat'];
                if exist(fullfile(vidDataDir, changeLogFile), 'file')
                    
                    load(fullfile(vidDataDir, changeLogFile)); % 'renameList', 'delList', 'newNames'
                    
                    % Extract ID numbers from filenames
                    disp('Extracting ID numbers...')
                    match_func = @(x) str2double(regexp(x, '(?<=tid_).*', 'match'));
                    vidDirNames = {trialVidDirs.name};
                    trialNums = cellfun(match_func, vidDirNames);
                    
                    % Fix sorting order to represent actual time of frame acquisition
                    disp('Sorting file names...')
                    [~, sortOrder] = sort(trialNums);
                    vidDirsSorted = vidDirNames(sortOrder);
                    
                    % Remove excess behavior video
                    disp('Moving behavior video...')
                    for iDir = 1:numel(trialVidDirs)
                        if iDir > numel(newNames)
                            movefile(fullfile(vidDataDir, vidDirsSorted{iDir}), fullfile(vidDataDir, ...
                                'BadTrialsBackup', vidDirsSorted{iDir}));
                        end
                    end
                    
                    % Double check that directories have actually been moved before continuing
                    loopCount = 0;
                    while exist(fullfile(vidDataDir, vidDirsSorted{iDir}), 'dir') && loopCount < 15
                        write_to_log('Waiting 5 seconds for move to complete...', mfilename);
                        write_to_log(fullfile(vidDataDir, vidDirsSorted{iDir}), mfilename);
                        loopCount = loopCount + 1;
                        pause(5);
                    end
                    
                    % Update stim metadata file
                    disp('Updating metadata file...')
                    mdFile = dir(fullfile(vidDataDir, ['metadata*sid_', num2str(sid), '*bid_', ...
                        num2str(bid), '.mat']));
                    load(fullfile(mdFile.folder, mdFile.name)); % --> 'metaData' with fields 'nTrials',
                    % 'trialDuration', 'stimTypes', 'sid',
                    % 'taskFile', 'outputData'
                    nTrialsOld = metaData.nTrials;
                    sampPerTrial = size(metaData.outputData, 1) / metaData.nTrials;
                    metaData.nTrials = numel(newNames);
                    metaData.stimTypes = metaData.stimTypes(1:numel(newNames));
                    metaData.outputData = metaData.outputData(1:sampPerTrial * metaData.nTrials, :);
                    movefile(fullfile(mdFile.folder, mdFile.name), fullfile(vidDataDir, 'BadTrialsBackup', ...
                        ['backup_', mdFile.name]));
                    save(fullfile(mdFile.folder, mdFile.name), 'metaData');
                    
                    % Update fictrac data files
                    disp('Updating FicTrac data file...')
                    ftFile = dir(fullfile(vidDataDir, ['fictracData*sid_', num2str(sid), '*bid_', ...
                        num2str(bid), '.mat']));
                    load(fullfile(ftFile.folder, ftFile.name)); % 'blockData' --> [nSamples, nChannels]
                    sampPerTrial = size(blockData, 1) / nTrialsOld;
                    blockData = blockData(1:sampPerTrial * numel(newNames), :);
                    movefile(fullfile(ftFile.folder, ftFile.name), fullfile(vidDataDir, ...
                        'BadTrialsBackup', ['backup_', ftFile.name]));
                    save(fullfile(ftFile.folder, ftFile.name), 'blockData');
                else
                    disp(['Change log file ', fullfile(vidDataDir, changeLogFile), ' not found!'])
                    write_to_log(['Change log file ', fullfile(vidDataDir, changeLogFile), ...
                        ' not found!'], mfilename)
                end%if
            else
                disp('Backup directory already contains files from this block!')
                write_to_log('Backup directory already contains files from this block!', mfilename);
            end%if
        end%iBlock
    end
catch ME
    write_to_log(getReport(ME), mfilename);
end%function