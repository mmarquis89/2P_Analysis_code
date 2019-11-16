
parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\20191112-1_38A11-Chrimson-6m\ProcessedData';

trialNum = 12;


% sz = size(raw);
% raw = raw(:, :, 1:(10 * floor(sz(3)/10)));
% dsRaw = squeeze(mean(reshape(raw, sz(1), sz(2), 10, []), 3));
% 
% sz = size(reg);
% reg = reg(:, :, 1:(10 * floor(sz(3)/10)));
% dsReg = squeeze(mean(reshape(reg, sz(1), sz(2), 10, []), 3));

refImgFiles = dir(fullfile(parentDir, 'refImages_reg_trial*.mat'));
allTrialRefImages = [];
for iFile = 1:19
   currFileName = refImgFiles(iFile).name;
   load(fullfile(parentDir, currFileName));
   allTrialRefImages(:,:, iFile) = refImages;
end

%%
lumThresh = 400;
for iTrial = 1:size(allTrialRefImages, 3)
   figure(iTrial); clf;
   imshow(allTrialRefImages(:,:,iTrial), [0 lumThresh]);   
end