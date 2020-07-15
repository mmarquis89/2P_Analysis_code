
%% Break data into chunks
testFl = a.analysisOutput.rawFlData;
% test = cell2padded_mat(currSubset.moveSpeed);
% test = test(:, 1) * 25 * 4.5;
% testSpeed = repeat_smooth(test, 15, 'dim', 1, 'smWin', 7);

testSpeed = abs(a.analysisOutput.fwSpeedData);
% disp(currSubset.originalTrialCount);

% chunkLength = 960;
nChunks = 3;
nBins = 15;
minX = 0;


% cm = [rgb('blue'); rgb('green'); rgb('red'); rgb('red'); rgb('green'); rgb('blue'); rgb('green'); ...
%         rgb('red'); rgb('red'); rgb('green');];
% cm = [rgb('green'); rgb('red'); rgb('red'); rgb('green'); rgb('green'); rgb('red'); ...
%         rgb('red'); rgb('green'); rgb('green');];
% cm = [rgb('green'); rgb('red'); rgb('green'); rgb('red')];
cm = [rgb('green'); rgb('red'); rgb('blue'); rgb('red'); rgb('red')];
% cm = repmat([rgb('red'); rgb('green')], 2, 1);


% nChunks = round((numel(testFl) / a.sourceDataTable.subset.volumeRate(1)) / chunkLength); 
chunkSize = floor(numel(testFl) / nChunks);
chunkedFl = reshape(testFl(1:(chunkSize * nChunks)), [], nChunks);
chunkedSpeed = reshape(testSpeed(1:(floor(numel(testSpeed) / nChunks)) * nChunks), [], nChunks);

figure(12);clf; imagesc(chunkedFl'); title('chunkedFl');
figure(13);clf; imagesc(chunkedSpeed'); title('chunkedSpeed');

figure(14);clf; imagesc(sort(chunkedFl)'); title('sorted chunkedFl');
figure(15);clf; imagesc(sort(chunkedSpeed)'); title('sorted chunkedSpeed');

% chunkedFl = chunkedFl(:, 1:end-1);
% chunkedSpeed = chunkedSpeed(:, 1:end - 1);

%  Calculate dF/F for each one independently 
chunkedFlSorted = sort(chunkedFl);
chunkBaselines = [];
for iChunk = 1:size(chunkedFl, 2)
    chunkBaselines(iChunk) = chunkedFlSorted(round(sum(~isnan(chunkedFl(:, iChunk))) * 0.05), iChunk);
end
base = repmat(chunkBaselines, size(chunkedFl, 1), 1);
chunkedDff = (chunkedFl - base) ./ chunkedFl;

% chunkedDff = chunkedFl;

figure(16);clf
imagesc(chunkedDff'); title('chunked dFF');
figure(17);clf;
histogram2(chunkedSpeed(:), chunkedDff(:), 50, 'displaystyle', 'tile')
figure(18); clf;hold on
if size(chunkedSpeed, 2) < 20
%     cm = jet(size(chunkedSpeed, 2)).*0.9;
    for i=1:size(chunkedDff, 2)
        [midpoints, means, SEM] = MoveSpeedAnalysis.bin_data(chunkedSpeed(:, i), chunkedDff(:, i), nBins, minX);
        errorbar(midpoints, means, SEM, '-o', 'color', cm(i, :));
    end
    xlabel('moveSpeed(mm/sec)')
    ylabel('mean binned dF/F')
end

% Calculate dF/F using a sliding window instead of static chunks
% base = mov_percentile(testFl, round(chunkSize / 2), 0.05)';
base = mov_percentile(testFl, 380, 0.05)';
dff = (testFl - base) ./ base;
figure(21);clf;
histogram2(testSpeed(~isnan(testSpeed)), dff(~isnan(dff)), 50, 'displaystyle', 'tile')
figure(22); clf;
[midpoints, means, SEM] = MoveSpeedAnalysis.bin_data(testSpeed, dff, nBins, minX);
errorbar(midpoints, means, SEM, '-o', 'color', 'k');
xlabel('moveSpeed')
ylabel('mean sliding dF/F')








