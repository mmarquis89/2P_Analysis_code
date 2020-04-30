
%% Blocks of processed imaging data

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017\2017_11_27_exp_2\sid_2';

imgDataFiles = dir(fullfile(parentDir, 'imagingData_reg_block*.mat'));
medVals = [];
minVals = [];
meanVals = [];
for iFile = 1:numel(imgDataFiles)
    disp(iFile)
    load(fullfile(parentDir, imgDataFiles(iFile).name), 'imgData'); % --> [y, x, plane, vol, trial]
    sz = size(imgData);
    imgDataPerm = permute(reshape(imgData, sz(1), sz(2), sz(3), []), [4 1 2 3]); % --> [expVol, y, x, plane]
    medVals = [medVals, median(reshape(imgDataPerm, size(imgDataPerm, 1), []), 2)'];
    minVals = [minVals, min(reshape(imgDataPerm, size(imgDataPerm, 1), []), [], 2)'];
    meanVals = [meanVals, mean(reshape(imgDataPerm, size(imgDataPerm, 1), []), 2)'];
end


%% Cdata .mat files

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\2017\2017_11_27_exp_2\sid_2';

planeNum = 6;

imgDataFiles = dir(fullfile(parentDir, 'cdata*.mat'));
medVals = [];
minVals = [];
meanVals = [];
for iFile = 1:numel(imgDataFiles)
    disp(iFile)
    load(fullfile(parentDir, imgDataFiles(iFile).name), 'tifData'); % --> [y, x, plane, vol]
    tifDataPerm = permute(tifData, [4 1 2 3]);                      % --> [vol, y, x, plane]
    medVals = [medVals, median(reshape(tifDataPerm, size(tifDataPerm, 1), []), 2)'];
    minVals = [minVals, min(reshape(tifDataPerm, size(tifDataPerm, 1), []), [], 2)'];
    meanVals = [meanVals, mean(reshape(tifDataPerm, size(tifDataPerm, 1), []), 2)'];
end

%%

sz = size(wholeSession); 
wholeSession = reshape(wholeSession, sz(1), sz(2), sz(3), []);  % --> [y, x, plane, expVol]
wholeSession = permute(wholeSession, [4 1 2 3]);                % --> [expVol, y, x, plane]
wholeSession = reshape(wholeSession, size(wholeSession, 1), []);% --> [expVol, px]
medSessionVals = median(wholeSession, 2);
minSessionVals = min(wholeSession, [], 2);
meanSessionVals = mean(wholeSession, 2);


%%

sz = size(wholeSession);
wholeSession = reshape(wholeSession, sz(1), sz(2), sz(3), []);  % --> [y, x, plane, expVol]
wholeSession = permute(wholeSession, [4 1 2 3]);                % --> [expVol, y, x, plane]
wholeSession = reshape(wholeSession, size(wholeSession, 1), []);% --> [expVol, px]
medRegVals = median(wholeSession, 2);
minRegVals = min(wholeSession, [], 2);
meanRegVals = mean(wholeSession, 2);

%%

xx = 1:129:size(meanVals, 2);

figure(1); clf; hold on
plot(minVals)
plot(medVals)
plot(smoothdata(meanVals, 'gaussian', 7))
plot(xx, ones(size(xx)) * 20, 'o', 'color', 'm')

% figure(2); clf; hold on
% plot(minSessionVals);
% plot(medSessionVals);
% plot(smoothdata(meanSessionVals, 'gaussian', 11));
% plot(xx, ones(size(xx)) * 50, 'o')
% 
% figure(3); clf; hold on
% plot(minRegVals);
% plot(medRegVals);
% plot(smoothdata(meanRegVals, 'gaussian', 11));
% plot(xx, ones(size(xx)) * 50, 'o')

%%

sz = size(imgData);
rsData = reshape(imgData, sz(1), sz(2), sz(3), []); % [y, x, plane, vol]
permData = permute(rsData, [4 1 2 3]); % [vol, y, x, plane]
sz = size(permData);
medPxByVol = median(reshape(permData, sz(1), []), 2); % [vol]
sz = size(rsData);
medPxRep = permute(repmat(medPxByVol, 1, sz(1), sz(2), sz(3)), [2 3 4 1]); % [y, x, plane, volume]

baseSubData = rsData - medPxRep;
meanVals = squeeze(multi_mean(baseSubData, [1 2]));        
sz = size(baseSubData);
medVals = squeeze(median(reshape(permute(baseSubData, [3 4 1 2]), sz(3), sz(4), []), 3)); % --> [plane, vol] 
minVals = squeeze(min(reshape(permute(baseSubData, [3 4 1 2]), sz(3), sz(4), []), [], 3));% --> [plane, vol] 
