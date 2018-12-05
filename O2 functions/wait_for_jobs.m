function [jobArr] = wait_for_jobs(jobArr)

% Pause MATLAB execution until all jobs in SLURM job cell array are complete
    for iTrial = 1:numel(jobArr)
        jobArr{iTrial}.wait;
    end
end