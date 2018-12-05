function pca_calc(parentDir, sessionDataFile)
%===================================================================================================
% 
% Run PCA on all the pixels from each plane and save the results as a separate .mat file
% 
%===================================================================================================
try
% Load session data
load(fullfile(parentDir, sessionDataFile)) % --> 'wholeSession'

sz = size(wholeSession);

% Run PCA on each plane
pcaData = zeros(sz(1), sz(2), sz(4)-1, sz(3)); 
explained = [];
for iPlane = 1:sz(3)
   currData = mean(squeeze(wholeSession(:,:,iPlane, :, :)),4); % --> [y, x, volume]
   disp(['pcaData: ', num2str(size(currData))])
   [n,m,d] = size(currData);
   data2D = double(reshape(currData, [(n*m), d]));               % --> [pixel, volume]
   disp(['data2D: ', num2str(size(data2D))])
   tData2D = data2D';                                           % --> [volume, pixel]
   disp(['tData2D: ', num2str(size(tData2D))])
   [coeff,~,~,~,ex] = pca(tData2D);
   disp(['coeff: ', num2str(size(coeff))])
   coeffReshaped = reshape(coeff, [n, m, d-1]);                 % --> [y, x, pc]
   disp(['coeffReshaped: ', num2str(size(coeffReshaped))])
   pcaData(:, :, :, iPlane) = coeffReshaped;                    % --> [y, x, pc, plane]
   explained(:, iPlane) = ex;                                   % --> [pc, iPlane] 
end

% Discard lower-value PCs if file size is large
if size(pcaData, 3) > 50
    pcaData = pcaData(:, :, 1:50, :);    
end

% Save output 
save(fullfile(parentDir, ['PCA_data_', sessionDataFile]), 'pcaData', 'explained', '-v7.3');
catch ME
    write_to_log(ME.message, mfilename);
end%try
end