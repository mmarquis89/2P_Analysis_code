function pcaTest()
addpath('/home/mjm60/HelperFunctions')
c = parcluster; 

parentDir = '/home/mjm60/'


% Load test data
write_to_log('Loading test data', mfilename);
load(fullfile(parentDir, 'testData.mat')); % --> sessionData
write_to_log('Test data loaded', mfilename);

% Create downsampled datasets for testing
dsFactors = [1];
dsData = []; testTimes = [];
for i = 1:numel(dsFactors)
    dsData{i} = imresize(sessionData, dsFactors(i));
    sz = size(dsData{i});
    dsData{i} = single(reshape(dsData{i}, [sz(1)*sz(2), sz(3)])');
%     sz = size(sessionData);
%     dsData{i} = single(reshape(sessionData, [sz(1)*sz(2), sz(3)])');

    % Run tests
    for iTest = 1:5
        write_to_log(['dsFactor: ', num2str(dsFactors(i)), '  test: ', num2str(iTest)], mfilename);
        tic
        [pcaData, ~,~,~,ex] = pca(dsData{i}(1:(1000*iTest), :));
        testTimes(i, iTest) = toc; % --> [imgSize, nVolumes]
    end
end

% Save output times
save(fullfile(parentDir, 'testTimes2.mat'), 'testTimes');

end