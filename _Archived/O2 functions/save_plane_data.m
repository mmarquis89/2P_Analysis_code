function save_plane_data(parentDir, sessionDataFile)

    addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

    % Load data
    [analysisMetadata, wholeSession] = load_imaging_data(parentDir, sessionDataFile);
    
    % Save one file containing the session data for each plane
    for iPlane = 1:analysisMetadata.nPlanes
        planeData = squeeze(wholeSession(:,:,iPlane,:,:)); % --> [y, x, volume, trial]
        fileName = ['plane_', num2str(iPlane), '_sessionData'];
        save(fullfile(parentDir, fileName), 'planeData', '-v7.3')
    end
end