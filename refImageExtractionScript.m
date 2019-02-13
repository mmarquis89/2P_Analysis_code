
parentDirs = { 'D:\Dropbox (HMS)\2P Data\Imaging Data\2018 Jan-Mar', ...
    'D:\Dropbox (HMS)\2P Data\Imaging Data\2018 Apr-May', ...
    'D:\Dropbox (HMS)\2P Data\Imaging Data', ...
    };

refImgData = [];

for iDir = 1:numel(parentDirs)
    parentDir = parentDirs{iDir};
    
    % Get list of all exp dirs in current parent dir
    expDirs = dir(fullfile(parentDir, '2018_*'));
    
    % Loop through each exp and each sid dir
    for iExp = 1:numel(expDirs)
        currExp = expDirs(iExp).name;
        disp(currExp);
        
        sidDirs = dir(fullfile(parentDir, currExp, 'sid_*'));
        for iSid = 1:numel(sidDirs)
            
            % Extract sid number
            currSidDir = fullfile(sidDirs(iSid).folder, sidDirs(iSid).name);
            sidNumCell = regexp(sidDirs(iSid).name, '(?<=sid_).', 'match');
            sid = str2double(sidNumCell{:});
            
            % Identify ref image files
            refImgFiles = dir(fullfile(currSidDir, '*refIm*.mat'));
            
            % Save reference image data for each file
            for iFile = 1:numel(refImgFiles)
                load(fullfile(refImgFiles(iFile).folder, refImgFiles(iFile).name), 'refImages');
                
                % Convert from cell to numeric array
                currRefImg = [];
                for iPlane = 1:numel(refImages)
                    currRefImg(:,:,iPlane) = refImages{iPlane};
                end
                
                refImgData(end + 1).expDate = currExp;
                refImgData(end).sid = sid;
                refImgData(end).fileName = refImgFiles(iFile).name;
                refImgData(end).refImages = currRefImg;
            end
        end
    end
end


%% Separate registered files

regFiles = zeros(numel(refImgData), 1);
for iFile = 1:numel(refImgData)
    
    if regexp(refImgData(iFile).fileName, '(Reg|sessionFile)') 
        regFiles(iFile) = 1;
    end
    
end

regRefImgData = refImgData(logical(regFiles));


%%

for iFile = [65 66]%1:numel(regRefImgData)
    
    currRefImg = regRefImgData(iFile).refImages;
    
    f = figure(iFile); clf;
    f.Position = [-1910 10 1900 980];
    for iPlane = 1:size(currRefImg, 3)
        subplot(3, 4, iPlane)
        minVal = min(as_vector(currRefImg(:,:,iPlane)));
        maxVal = max(as_vector(currRefImg(:,:,iPlane)));
        dispRange = [minVal, maxVal * 0.8];
        if dispRange(2) <= dispRange(1)
           dispRange(2) = maxVal; 
        end
        imshow(currRefImg(:,:,iPlane), dispRange);
    end
    
    figTitle = regexprep([regRefImgData(iFile).expDate, '_sid_', num2str(regRefImgData(iFile).sid)], '\_', '\\_');
    suptitle(figTitle);
    
    
end

%%

goodCandidates = {'2018_03_23_exp_2', '2018_10_09_exp_2', '2018_10_18_exp_2', '2018_10_24_exp_3', ...
    '2018_10_30_exp_1', '2018_11_11_exp_3'};

secondaryCandidates = {'2018_01_27_exp_2', '2018_02_07_exp_2', '2018_07_07_exp_1', '2018_10_24_exp_2', ...
    '2018_10_30_exp_2', '2018_11_07_exp_1', '2018_11_07_exp_2', '2018_11_09_exp_1', '2018_11_09_exp_2', ...
    '2018_11_11_exp_1', '2018_11_20_exp_1'};

goodCandidateInds = ismember({regRefImgData.expDate}, goodCandidates);
goodCandidateRefImg = {regRefImgData(goodCandidateInds).refImages};

secCandidateInds = ismember({regRefImgData.expDate}, secondaryCandidates);
secCandidateRefImg = {regRefImgData(secCandidateInds).refImages};




