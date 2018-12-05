

parentDir = 'B:\Dropbox (HMS)\2P Data\Behavior Vids\2018_05_25\_Movies\FicTracData';
dirFiles = dir(fullfile(parentDir, 'sid*.dat'));

rowCounts = []; frameCounts = []; ftFrameCounts = [];
for iFile = 1:length(dirFiles)
test = csvread(fullfile(parentDir, dirFiles(iFile).name));
rowCounts(iFile,1) = size(test,1);
frameCounts(iFile,1) = test(end, 1);
ftFrameCounts(iFile,1) = test(end, 23);
end
comp = [rowCounts, frameCounts, ftFrameCounts];
figure(1); imagesc(comp)