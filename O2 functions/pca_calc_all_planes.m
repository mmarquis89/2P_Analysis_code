function pca_calc_all_planes(parentDir, sessionDataFile)
%===================================================================================================
% 
% Run PCA on all the pixels from each plane and save the results as a separate .mat file
% 
%===================================================================================================

% Load session data
load(fullfile(parentDir, sessionDataFile)) % --> 'wholeSession'

%% Run PCA on entire dataset
sz = size(wholeSession);
allVolSession = reshape(wholeSession, [sz(1), sz(2), sz(3), (sz(4)*sz(5))]);    % --> [y, x, plane, trialVolume]
sz = size(allVolSession);
data2D = double(reshape(allVolSession, [(sz(1)*sz(2)*sz(3)), sz(4)]));          % --> [pixel, trialVolume]
tData2D = data2D';                                                              % --> [trialVolume, pixel]
[coeff,score,latent,~,explained] = pca(tData2D);                                % [pixel, pc]
pcaData = reshape(coeff, [sz(1:3), sz(4)-1]);                                   % [y, x, plane, pc]

%%


% Save output 
save(fullfile(parentDir, ['PCA_data_all_planes', sessionDataFile]), 'pcaData', 'score', 'latent', 'explained', '-v7.3');

end%function

% 
% singlePlaneData = squeeze(wholeSession(:,:,6,:,:));
% 
% %%
% testPlane = 6;
% 
% 
% 
% 
% sz = size(singlePlaneData);
% allVolData = reshape(singlePlaneData, [sz(1:2), sz(3)*sz(4)]);  % --> [y, x, trialVolume]
% sz = size(allVolData);
% data2D = double(reshape(allVolData, [(sz(1)*sz(2)), sz(3)]));   % --> [pixel, trialVolume]
% tData2D = data2D';                                              % --> [trialVolume, pixel]
% % [coeff2,latent2,explained2] = pcacov(tData2D);                   % --> [pixel, pc]
% % [coeff,score,latent,tSquared,explained] = pca(tData2D);        % --> [pixel, pc]
% % pcaData = reshape(coeff, [sz(1:2), sz(3)-1]);                   % --> [y, x, pc]
%   pcaDataCov = reshape(coeff2, [sz(1:2), sz(3)-1]);                   % --> [y, x, pc]
%   
% figure(1);clf;
% plot(explained2(1:100))
% 
% %%
% testPC = 28;
% 
% figure(2);clf
% imagesc(pcaData(:,:,testPC));
% colormap(gca, 'bluewhitered');
% 
% 
% %% Discard uninformative PCs 
% 
% 
% %% Recreate original data with limited PCs
% mu = mean(tData2D);
% sigma = std(tData2D);
% nPCs = [1:10];
% Xhat = score(:, nPCs) * coeff(:, nPCs)';
% 
% % Xhat = bsxfun(@times, Xhat, sigma);
% Xhat = bsxfun(@plus, Xhat, mu);
% 
% %%
% 
% newData = Xhat';
% sz = size(newData);
% newData = reshape(newData, [128, 256, sz(2)]);
% implay(movmean(newData, 3, 3))
% % newRefImg = mean(newData, 3);
% % figure; imshow(newRefImg, [0 800])
% 
% %%
% 
% trimData = coeff(:,1:5);
% 
% q = 3;
% Mdl = rica(trimData, q, 'VerbosityLevel', 1);
% 
% newX = transform(Mdl,trimData); 






