function normcorre_registration(parentDir, fileName, varargin)
%===================================================================================================
%
% Uses the NoRMCorre algorithm to apply a rigid motion correction to pre-processed imaging data. It 
% must be saved as a .mat file containing the variables 'trialType', 'origFileNames', 'expDate', and 
% either 'wholeSession' or 'regProduct' (the latter in case I want to register files that have 
% already been processed by my old registration algorithm). The imaging data array (which must have 
% dimensions [y, x, plane, volume, trial])is the only variable that is actually used, the rest are 
% just re-saved as a new .mat file in the same directory as the source file.
% 
%  INPUT:
%       parentDir   = the full path to the directory containing the data to be registered
%       fileName    = the name of the input .mat file (this will also affect the output file name)
%
% OPTIONAL NAME-VALUE PAIR INPUTS:
%       
%       'OutputDir' = (default: parentDir) the directory to save the processed data in
%===================================================================================================

try
write_to_log('Starting NoRMCorre registration...', mfilename)

addpath('/home/mjm60/HelperFunctions') % if running on O2 cluster
addpath('/home/mjm60/NoRMCorre-master/NoRMCorre-master') % if running on O2 cluster

% Parse optional arguments
p = inputParser;
addParameter(p, 'OutputDir', parentDir);
parse(p, varargin{:});
outputDir = p.Results.OutputDir;

write_to_log('Loading array...', mfilename);

% Load data
load(fullfile(parentDir, fileName)); % 'sessionData' % --> [x, y, volume]

write_to_log('Session loaded', mfilename);

% Make session folder for new files if necessary
if ~isdir(outputDir)
    mkdir(outputDir)
end

% Set parameters
options_rigid = NoRMCorreSetParms('d1', size(sessionData, 1), 'd2', size(sessionData, 2), ...
                    'max_shift', [25, 25], ...
                    'init_batch', 100, ...
                    'us_fac', 50 ...
                    ); 
% options_nonRigid = NoRMCorreSetParms('d1', size(concatSession, 1), 'd2', size(concatSession, 2), 'd3', size(concatSession, 3), ...
%                     'max_shift', [25, 25, 2], ...
%                     'init_batch', 100, ...
%                     'grid_size', [100, 100] ...
%                     );
                
write_to_log('Starting registration...', mfilename);

% Rigid registration
tic; [sessionData, ~, regTemplate, ~] = normcorre_batch(sessionData, options_rigid); toc

% Save registered data
save(fullfile(parentDir, ['rigid_reg_', fileName]), ...
        'sessionData', '-v7.3')

% Save reference images
refImg = mean(sessionData, 3);
sz = size(sessionData);
sessionData = reshape(sessionData, sz(1), sz(2), 100, []);
timePointRefImg = squeeze(mean(sessionData, 3)); % --> [y, x, section of volumes]
save(fullfile(parentDir, ['refImage_reg_', fileName]), 'refImg', ...
        'regTemplate', 'timePointRefImg', '-v7.3');

catch ME
    write_to_log(getReport(ME), mfilename);
end%try
end