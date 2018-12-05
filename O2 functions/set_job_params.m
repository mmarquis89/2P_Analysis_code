function c = set_job_params(c, queueName, timeLimit, memGB, jobName)
%===================================================================================================
% Updates resource request parameters for SLURM jobs
%
% INPUTS:
%
%       c         = the SLURM cluster object that you are using to run your jobs  
% 
%       queueName = the SLURM queue to submit to ('priority', 'short', 'medium', 'long')
% 
%       timeLimit = time limit for the jobs in minutes
% 
%       memGB     = the memory in GB to be requested
%
%       jobName   = the name of the job in SLURM
% 
%===================================================================================================

    c.AdditionalProperties.WallTime = num2str(timeLimit);
    c.AdditionalProperties.QueueName = queueName;
    c.AdditionalProperties.AdditionalSubmitArgs = ['-c 1 --mem=', num2str(memGB), 'G --job-name=', jobName];
    
end