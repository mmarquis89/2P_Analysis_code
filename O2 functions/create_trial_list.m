function create_trial_list(parentDir, sid)

% Create a comma-separated text file listing all the folders to be processed 
% for a specific session and another one listing the tid of each dir

disp(parentDir)

dirContents = dir(fullfile(parentDir, '*sid*tid*'));
dirFile = fopen(fullfile(parentDir, ['sid_', num2str(sid), '_dirList.txt']), 'w');
disp(num2str(numel(dirContents)))
tid = 0;
for i = 1:numel(dirContents)
    disp(fullfile(parentDir, dirContents(i).name))
    if isdir(fullfile(parentDir, dirContents(i).name))
        currSid = str2double(regexp(dirContents(i).name, '(?<=sid_).*(?=_tid)', 'match'));
        if currSid == sid
            tid = tid + 1;
            writeStr = [fullfile(parentDir, [dirContents(i).name]), ',', num2str(tid), '\n'];
            fprintf(dirFile, writeStr);
        end
    end
end
fclose(dirFile);

