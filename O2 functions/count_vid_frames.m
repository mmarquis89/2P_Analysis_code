function frameCounts = count_vid_frames(parentDir, sid)
try
addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster

logFile = fopen(fullfile(parentDir, ['sid_', num2str(sid), '_frameCounts.txt']), 'w');
vidFiles = dir(fullfile(parentDir, ['sid_', num2str(sid), '_tid*.avi']));

frameCounts = [];
for iFile = 1:numel(vidFiles)
    frameCounts(iFile) = count_frames(fullfile(parentDir, vidFiles(iFile).name));
    fprintf(logFile, [num2str(frameCounts(iFile)), ',', vidFiles(iFile).name, '\n']);
end
fclose(logFile);
catch ME
    write_to_log(ME.message, mfilename);
end%try
end