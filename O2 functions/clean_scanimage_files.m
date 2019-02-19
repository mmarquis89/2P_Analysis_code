function clean_scanimage_files(imgDataDir, sid)


try
%     N_PLANES = 16; % NOTE HARDCODED VALUE!
    
    if ~exist(fullfile(imgDataDir, ['sid_', num2str(sid), '_volume_counts.mat']), 'file')
        
        % Count and save number of volumes in imaging data files
        disp('Counting files...')
        trialFiles = dir(fullfile(imgDataDir, ['*sid_', num2str(sid), '*.tif']));
        bids = str2double(regexp({trialFiles.name}, '(?<=bid_).*(?=_dur)', 'match', 'once'));
        bidList = unique(bids);
        volCounts = []; frameCounts = [];
        for iTrial = 1:numel(trialFiles)
%             write_to_log(trialFiles(iTrial).name, mfilename)
            [currData, tifMetadata] = read_patterned_tifdata(fullfile(imgDataDir, ...
                        trialFiles(iTrial).name)); % --> [y, x, allVols]
%             currData = permute(currData, [2 1 3]);
            sz = size(currData);
%             write_to_log(num2str(sz), mfilename);
            nPlanes = frameStringKeyLookup(tifMetadata.tifinfo.Software, 'hStackManager.numSlices');
            nFlybackFrames = frameStringKeyLookup(tifMetadata.tifinfo.Software, ...
                        'hFastZ.numDiscardFlybackFrames');
            frameCounts(iTrial) = sz(3);
            volCounts(iTrial) = sz(3) / (nPlanes + nFlybackFrames);
%             write_to_log(['VolCounts = ', num2str(volCounts(iTrial))], mfilename);
        end
        disp(['Saving ', fullfile(imgDataDir, ['sid_', num2str(sid), '_volume_counts.mat'])])
        
        save(fullfile(imgDataDir, ['sid_', num2str(sid), '_volume_counts.mat']), 'volCounts', ...
            'trialFiles', 'bidList', '-v7.3');
        
        
        % Fix files
        disp('Fixing file names...')
        for iBlock = 1:numel(bidList)
            currVolCounts = volCounts(bids == bidList(iBlock));
            currFiles = trialFiles(bids == bidList(iBlock));
            
            
            % Check whether there are single-volume trials
            if min(currVolCounts) == 1
                
                disp('Eliminating single-volume trials...')
                write_to_log('Eliminating single-volume trials...');
                
                targetVols = mode(currVolCounts);
                renameList = {currFiles(~(currVolCounts < targetVols)).name}';
                delList = {currFiles(currVolCounts < targetVols).name}';
                newNames = {currFiles(1:numel(renameList)).name}';
                if ~isdir(fullfile(imgDataDir, 'BadTrialsBackup'))
                    mkdir(fullfile(imgDataDir, 'BadTrialsBackup'))
                end
                
                % Save record of changes
                save(fullfile(imgDataDir, ['sid_', num2str(sid), '_bid_', num2str(bidList(iBlock)), ...
                    '_single_vol_trials.mat']), 'renameList', 'delList', 'newNames')
                
                % Remove single-volume trials
                for iFile = 1:numel(delList)
                    try
                        movefile(fullfile(imgDataDir, delList{iFile}), fullfile(imgDataDir, ...
                                'BadTrialsBackup', delList{iFile}));
                    catch
                        pause(1);
                        movefile(fullfile(imgDataDir, delList{iFile}), fullfile(imgDataDir, ...
                                'BadTrialsBackup', delList{iFile}));
                    end
                end
                
                % Rename other trials as needed
                for iFile = 1:numel(renameList)
                    oldName = fullfile(imgDataDir, renameList{iFile});
                    newName = fullfile(imgDataDir, newNames{iFile});
                    if ~strcmp(oldName, newName)
                        try
                            movefile(oldName, newName);
                        catch
                            pause(1);
                            movefile(oldName, newName);
                        end
                    end
                end%for
            end%if
            
            
            % Check whether the total number of volumes fluctuates
            modeVols = mode(currVolCounts);
            remVols = mode(currVolCounts(currVolCounts ~=modeVols));
            if abs(modeVols - remVols) < 2
                
                disp('Trimming extra volumes off ends of trials...')
                write_to_log('Trimming extra volumes off ends of trials...', mfilename);
                
                % Trim one volume off all of the longer trials and re-save as .mat file
                minVols = min([modeVols, remVols]);
                for iTrial = 1:numel(currFiles)
                    currFile = fullfile(imgDataDir, currFiles(iTrial).name);
                    disp(currFile)
                    % currData = permute(read_patterned_tifdata(currFile), [2 1 3]);
                    currData = read_patterned_tifdata(currFile);
                    sz = size(currData);
                    nVolumes = sz(3) / (nPlanes + nFlybackFrames);
                    rsData = reshape(currData, sz(1), sz(2), (nPlanes + nFlybackFrames), nVolumes); % --> [y, x, plane, volume]
                    tifData = rsData(:,:,:, 1:minVols);
                    saveName = regexprep(currFile, '.tif', '.mat');
                    save(saveName, 'tifData', '-v7.3');
                end
            end
            
        end%iBlock
    else
        write_to_log([fullfile(imgDataDir, ['sid_', num2str(sid), '_volume_counts.mat']), ...
                ' already exists!'], mfilename);
    end%if
catch ME
    disp(getReport(ME))
    write_to_log(getReport(ME), mfilename);
end%function