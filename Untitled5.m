
roiDataCat = [];

% Create matrix of raw data
for iROI = 1:numel(currROIData)
   roiDataCat(iROI, :) = currROIData(iROI).data; 
end


% Calculate dF/F
roiDataSort = sort(roiDataCat, 2);
baselineF = median(roiDataSort(:, 1:(size(roiDataSort, 2) * 0.05)), 2);
baselineMat = repmat(baselineF, 1, size(roiDataCat, 2));
roiDff = (roiDataCat - baselineMat) ./ baselineMat;
roiDff = roiDff'; % --> [volume, ROI]

% Get z-score
roiZScore = zscore(roiDff); % --> [volume, ROI]


% Calculate power spectrum and phase
transform = fft(roiDff, 16, 1);
transform(1,:) = [];
n = 16;
phase = angle(transform(1:8, :));
phaseValue = squeeze(phase(2, :));

s = .5.^(1:0.25:5);
pxx = periodogram(roiDff, [], s, 1, 'power');
posEight = find(s == .125);
powerValue = pxx(posEight, :);


% Combine dF/F from L&R glomeruli and calculate PVA
left_dff = roiDff(1:8,:);
right_dff = roiDff(9:16,:);
mean_dff = (left_dff + right_dff)./2;
mu = circ_mean(repmat([-7*pi/8:pi/4:7*pi/8], size(mean_dff,2),1), mean_dff', 2);
dff_pva = ((mu)/pi*4+4.5);
dff_pva = dff_pva';
dff_pva_rad = mu';

%%

figure(1);clf;imagesc(roiDff(:, 1:8)'+ roiDff(:, 9:end)');

%%
pq = (numel(panelsPos) / nVolumes);
nVolumes = size(roiDff, 1);
x = 1:numel(panelsPos);
v = panelsPos;
xq = pq:pq:numel(panelsPos); 
test = interp1(x, v, xq, 'nearest');



figure(10);clf;plot(panelsPos(round(xq)));




