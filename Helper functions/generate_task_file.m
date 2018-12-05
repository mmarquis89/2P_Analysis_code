function generate_task_file(tasks, reps, savePath, fileName)
%===================================================================================================
% 
% Generates a pseudo-randomized task file for 2P experiments, allowing user to customize number of 
% tasks and number of appearances of each task. The output list is saved as a .txt file, overwriting 
% any previous file with that name.
%
% INPUTS:
%       tasks     =  Nx1 cell array containing the name of one task in each cell.
%       
%       reps      =  Nx1 vector containing the number of times each task should appear in the file.
%
%       savePath  =  Full path to the directory where the task file is to be saved. Can pass [] to 
%                    generate a pop-up window to select a directory.
%
%       fileName  =  The task file name (without the .txt extension). Pass [] to use the default 
%                    naming convention of each task type and number in order concatenated with 
%                    underscores (e.g. "taskA_3_taskB_5").
%
%===================================================================================================

% Prompt user for save directory if none was provided
if isempty(savePath)
   savePath = uigetdir('C:\Users\wilson_lab\Desktop\Michael', 'Select a save directory'); 
end

% Use default file naming system if no name was provided
if isempty(fileName)
    nameStr = '';
    for iTask = 1:length(tasks)
        nameStr = [nameStr, tasks{iTask}, '_', num2str(reps(iTask)), '_'];
    end
    fileName = nameStr(1:end-1);
end

% Generate list of tasks
taskRep = [];
for iTask = 1:length(tasks)
   taskRep = [taskRep; repmat(tasks(iTask), reps(iTask), 1)];
end

% Randomize task list
taskRepRand = taskRep(randperm(length(taskRep)));

% Write list to file
taskFile = fopen(fullfile(savePath, [fileName, '.txt']), 'w+');
for iTask = 1:length(taskRepRand)
    fprintf(taskFile, '%s\r\n', taskRepRand{iTask});
end
fclose(taskFile);

end