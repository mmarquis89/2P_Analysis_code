function single_plane_pca_calc(parentDir, sessionDataFile)
%===================================================================================================
% 
% Run PCA on all the pixels from a single plane and save the results as a separate .mat file
% 
% 'sessionDataFile' should be an array with dimensions representing: [y, x, volume]
% 
%===================================================================================================
try
% Load session data
load(fullfile(parentDir, sessionDataFile)) % --> 'sessionData'

% Run PCA
ex = [];
[n, m, d] = size(sessionData);
sessionData = reshape(sessionData, [(n*m), d]); % --> [pixel, volume]
sessionData = sessionData';                     % --> [volume, pixel]
sessionData = single(sessionData);              % --> [volume, pixel]           
[coeffs,~,~,~, ex] = pca(sessionData);
coeffs = reshape(coeffs, [n, m, d-1]);          % --> [y, x, pc]

% Discard lower-ranked PCs if file size is large
if size(coeffs, 3) > 15
    coeffs = coeffs(:, :, 1:15, :);    
end

% Save output 
save(fullfile(parentDir, ['PCA_', sessionDataFile]), 'coeffs', 'ex', '-v7.3');

catch ME
    write_to_log(getReport(ME), mfilename);
end%try
end