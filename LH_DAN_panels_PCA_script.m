

startDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data';
parentDir = uigetdir(startDir, 'Select an experiment directory');


%% Run PCA on entire dataset


trialNum = 7;

load(fullfile(parentDir, ['imagingData_reg_trial_', ...
        pad(num2str(trialNum), 3, 'left', '0'), '.mat']), 'imgData')

sz = size(imgData);
data2D = double(reshape(imgData, [(sz(1)*sz(2)*sz(3)), sz(4)]));     % --> [pixel, trialVolume]
tData2D = data2D';                                                   % --> [trialVolume, pixel]
tic
[coeff,score,latent,~,explained] = pca(tData2D);                     % [pixel, pc]
tc = toc
pcaData = reshape(coeff, [sz(1:3), sz(4)-1]);                        % [y, x, plane, pc]

pcaData = pcaData(:, :, :, 1:100);
fileName = ['pcaData_trial_', pad(num2str(trialNum), 3, 'left', '0'), '.mat'];
save(fullfile(parentDir, fileName), 'pcaData')

%%

trialNum = 6;

load(fullfile(parentDir, ['pcaData_trial_', pad(num2str(trialNum), 3, 'left', '0'), '.mat']), ...
        'pcaData')
    
load(fullfile(parentDir, ['refImages_reg_trial_', ...
        pad(num2str(trialNum), 3, 'left', '0'), '.mat']), 'refImages')
    
%% 

nPCs = 7;
offset = 0;
   
nPlanes = size(pcaData, 3);
f = figure(1); clf;
f.Color = [1 1 1];
subaxis(nPCs + 1, nPlanes, 1, 1, ...
        'm', 0.01, ...
        'p', 0, ...
        'sh', 0.001, 'sv', 0)
for iPlane = 1:nPlanes
    
    % Plot reference images
    subaxis(nPCs + 1, nPlanes, iPlane, 1)
    imshow(refImages(:, :, iPlane), [0 700]);
    colormap(gca, 'gray');
    
    % Plot top 5 PCs
    for iPC = 1:nPCs
       subaxis(nPCs + 1, nPlanes, iPlane, iPC + 1) 
       imagesc(imgaussfilt(squeeze(pcaData(:, :, iPlane, iPC + offset)), 0.5)); 
       colormap(gca, 'bluewhitered'); 
       axis image
       axis off
    end
    

end


%%
figure(2); clf;
subplot(2,2,1)
imshow(refImages(:, :, planeNum),[0 700])
colormap(gca, 'gray')
% colormap('parula')
for iPlot = 2:4
    subplot(2, 2, iPlot); imagesc(imgaussfilt(planeData(:,:,iPlot-1),0.5));
    colormap(gca, 'bluewhitered')
end

figure(3); clf;
subplot(2,2,1)
for iPlot = 1:4
    subplot(2, 2, iPlot); imagesc(imgaussfilt(planeData(:,:,iPlot+3)));
    colormap(gca, 'bluewhitered')
end