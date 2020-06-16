
%% Break data into chunks
testFl = a.analysisOutput.flData;
% test = cell2padded_mat(currSubset.moveSpeed);
% test = test(:, 1) * 25 * 4.5;
% testSpeed = repeat_smooth(test, 15, 'dim', 1, 'smWin', 7);

testSpeed = a.analysisOutput.speedData;
% disp(currSubset.originalTrialCount);


% nChunks = 10;
nBins = 30;
minX = 0;

nChunks = round((numel(testFl) / a.sourceDataTable.subset.volumeRate(1)) / 30); 
chunkSize = floor(numel(testFl) / nChunks);
chunkedFl = reshape(testFl(1:(chunkSize * nChunks)), [], nChunks);
chunkedSpeed = reshape(testSpeed(1:(floor(numel(testSpeed) / nChunks)) * nChunks), [], nChunks);

figure(12);clf; imagesc(chunkedFl');
figure(13);clf; imagesc(chunkedSpeed');

figure(14);clf; imagesc(sort(chunkedFl)')
figure(15);clf; imagesc(sort(chunkedSpeed)');

%  Calculate dF/F for each one independently 
chunkedFlSorted = sort(chunkedFl);
chunkBaselines = chunkedFlSorted(round(size(chunkedFl, 1) * 0.05), :);
base = repmat(chunkBaselines, size(chunkedFl, 1), 1);
chunkedDff = (chunkedFl - base) ./ chunkedFl;
figure(16);clf
imagesc(chunkedDff')
figure(17);clf;
histogram2(chunkedSpeed(:), chunkedDff(:), 50, 'displaystyle', 'tile')
figure(18); clf;
[midpoints, means, SEM] = bin_data(chunkedSpeed(:), chunkedDff(:), nBins, minX);
errorbar(midpoints, means, SEM, '-o', 'color', 'k');
xlabel('moveSpeed(mm/sec)')
ylabel('mean binned dF/F')

% Calculate dF/F using a sliding window instead of static chunks
base = mov_percentile(testFl, round(chunkSize / 2), 0.05)';
dff = (testFl - base) ./ base;
figure(21);clf;
histogram2(testSpeed, dff, 50, 'displaystyle', 'tile')
figure(22); clf;
[midpoints, means, SEM] = bin_data(testSpeed, dff, nBins, minX);
errorbar(midpoints, means, SEM, '-o', 'color', 'k');
xlabel('moveSpeed')
ylabel('mean sliding dF/F')
% 
% Normalize both fl and speed data for each chunk
chunkedFlNorm = chunkedFl ./ (repmat(max(chunkedFl), size(chunkedFl, 1), 1));
chunkedSpeedNorm = chunkedSpeed ./ (repmat(max(chunkedSpeed), size(chunkedSpeed, 1), 1));

figure(19);clf; imagesc(chunkedFlNorm');
figure(20);clf; imagesc(chunkedSpeedNorm');

figure(21);clf;
histogram2(chunkedSpeedNorm(:), chunkedFlNorm(:), 50, 'displaystyle', 'tile')
figure(22); clf;
[midpoints, means, SEM] = bin_data(chunkedSpeedNorm(:), chunkedFlNorm(:), nBins, minX);
errorbar(midpoints, means, SEM, '-o', 'color', 'k');
xlabel('normalized moveSpeed')
ylabel('mean binned normalized F')






