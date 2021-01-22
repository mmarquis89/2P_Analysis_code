function params = load_plotting_params(expDir, varargin)
% Loads saved plotting parameters for the current experiment, if they exist.
% Optionally can be passed a name-value argument 'FileName' to override the default
% naming of 'plotting_params.mat'.
%
% Plotting parameters include:
%       stimNames
%       stimTrialGroups
%       stimGroupNames
%       stimEpochs
%       stimShadingColors
%       stimEpochNames
%       groupBounds
%       blockNames
%       blockShading

p = inputParser;
addParameter(p, 'FileName', 'plotting_params.mat');
parse(p, varargin{:});
fileName = p.Results.FileName;

if exist(fullfile(expDir, fileName), 'file')
    params = load(fullfile(expDir, fileName)); 
    disp('Loading plotting params...')
    disp('--------------------------------------------------------------------')
    disp(['stimNames = {', sprintf(' %s', params.stimNames{:}), '}'])
    disp('-----')
    disp(['stimTrialGroups = [', num2str(params.stimTrialGroups), ']'])
    disp('-----')
    disp(['stimGroupNames = {', sprintf(' %s', params.stimGroupNames{:}), '}'])
    disp('-----')
    disp('stimEpochs = '); disp(params.stimEpochs);
    disp('-----')
    disp('stimShadingColors = '), disp(params.stimShadingColors)
    disp('-----')
    disp('stimEpochNames = '); disp(params.stimEpochNames);
    disp('-----')
    disp(['groupBounds = ', num2str(params.groupBounds)])
    disp('-----')
    disp('blocknames = '), disp(params.blockNames);
    disp('-----')
    disp('blockShading = '), disp(params.blockShading);
    disp('--------------------------------------------------------------------')
else
    disp('No plotting parameter file found!')
    params = struct();
end

end