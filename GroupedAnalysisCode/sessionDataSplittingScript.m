

parentDir = 'F:\ImagingData';
expDirs = dir(fullfile(parentDir, '2018_02_09_exp_*'));
expDirs = expDirs([expDirs.isdir]);
expDirs = expDirs([1]); % 

for iExp = 1:numel(expDirs)
    sidDir = dir(fullfile(expDirs(iExp).folder, expDirs(iExp).name, 's*'));
    for iSid = 2%:numel(sidDir)
%         sessionDataFile = dir(fullfile(sidDir(1).folder, sidDir(1).name, 'rigid*'));
        sessionDataFile = dir(fullfile(sidDir(iSid).folder, sidDir(iSid).name, 'rigid*sessionFile*.mat'));
        if numel(sessionDataFile) > 0
            m = matfile(fullfile(sessionDataFile(1).folder, sessionDataFile(1).name));
            if contains('regProduct', fieldnames(m))
                sz = size(m, 'regProduct');
                oldName = 1;
            else
                sz = size(m, 'wholeSession');
                oldName = 0;
            end
            if sz(5) <= 80
                tic
                disp('Loading entire session file...')
                if oldName == 0
                    load(fullfile(sessionDataFile(1).folder, sessionDataFile(1).name), 'wholeSession')
                else
                    load(fullfile(sessionDataFile(1).folder, sessionDataFile(1).name), 'regProduct')
                end
                disp(['Entire session loaded in ', num2str(toc, 3), ' sec'])
            end
            
            blockStartTrials = 1:20:sz(5);
            nBlocks = numel(blockStartTrials);
            for iBlock = 1:nBlocks
                disp([expDirs(iExp).name, ', block ', num2str(iBlock), ' of ', num2str(nBlocks)])
                if iBlock < nBlocks
                    currBlockTrials = blockStartTrials(iBlock):(blockStartTrials(iBlock + 1) - 1);
                else
                    currBlockTrials = blockStartTrials(iBlock):sz(5);
                end
                if sz(5) <= 80
                    if oldName == 0
                        imgData = wholeSession(:, :, :, :, currBlockTrials);
                    else
                        imgData = regProduct(:, :, :, :, currBlockTrials);
                    end
                else
                    tic
                    disp('Loading single block from matfile...')
                    if oldName == 0
                        imgData = m.wholeSession(1:sz(1), 1:sz(2), 1:sz(3), 1:sz(4), currBlockTrials);
                    else
                        imgData = m.regProduct(1:sz(1), 1:sz(2), 1:sz(3), 1:sz(4), currBlockTrials);
                    end
                    disp(['Block loaded in ', num2str(toc, 3), ' sec'])
                end
                refImages = multi_mean(imgData, [4, 5]);
                
                tic
                save(fullfile(sessionDataFile(1).folder, ['imagingData_reg_block_', pad(num2str(iBlock), 2, ...
                    'left', '0'), '.mat']), 'imgData');
                disp(['Block saved to disk in ', num2str(toc, 3), ' sec'])
                
                save(fullfile(sessionDataFile(1).folder, ['refImages_reg_block_', pad(num2str(iBlock), 2, ...
                    'left', '0'), '.mat']), 'refImages');
                
                clear imgData
                
            end%iBlock
            
            clear wholeSession m regProduct
            
        end
    end
end%iExp