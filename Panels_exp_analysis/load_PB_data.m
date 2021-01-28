function [expMd, trialMd, roiData, ftData, flailingEvents, locEvents, panelsMetadata, wedgeData, ...
        glomData] = load_PB_data(parentDir, varargin)
% ==================================================================================================   
%  Loads all PB imaging data from multiple experiments from a grouped analysis diectory (can 
%  optionally pass a list of expIDs to use instead of loading the entire contents of the directory).
%  
%  INPUTS: 
%       parentDir           = full path to the grouped analysis directory full of data files.
%
%       expList (optional)  = cell array of a subset of expIDs to load from the directory. 
%
%
%  OUTPUTS:
%           
%       expMd               = experiment metadata table
% 
%       trialMd             = trial metadata table
% 
%       roiData             = ROI data table
% 
%       ftData              = FicTrac data table
% 
%       flailingEvents      = flailing event object
%
%       locEvents           = locomotion event object
% 
%       panelsMetadata      = panels metadata table
% 
%       wedgeData           = table containing the average fluorescence data and population vector 
%                             average (PVA) data for each "EB wedge" (i.e. pair of matching PB ROIs)
% 
%       glomData            = similar to wedgeData, but contains florescence data for each 
%                             individual glomerulus ROI instead of grouping into EB wedges
%
% ==================================================================================================

% Parse/process input arguments
if nargin < 2 || isempty(varargin{1})
    fileList = dir(fullfile(parentDir, '20*'));
    expList = unique(cellfun(@(x) x(1:10), {fileList.name}, 'uniformoutput', 0));
else
    expList = varargin{1};
end
    
% Load metadata 
[expMd, trialMd] = load_metadata(expList, parentDir);

% Load imaging data
roiData = load_roi_data(expList, parentDir);

% Load FicTrac data
ftData = load_ft_data(expList, parentDir);

% Load flailing event data (if any exists)
try
    eventData = load_event_data(expList, parentDir);
    flailingEvents = eventData.flailing;
catch
    flailingEvents = [];
end

% Load locomotion event data (if any exists)
try
    eventData = load_event_data(expList, parentDir);
    locEvents = eventData.locomotion;
catch
    locEvents = [];
end

% Load panels metadata
panelsMetadata = load_panels_metadata(expList, parentDir);

% ---------- Consolidate wedge data and calculate bump population vector averages -----------
    
% Calc dF/F for EB wedge ROIs
currData = roiData(~cellfun(@isempty, regexp(roiData.roiName, 'EB-', 'match')), :);
currData.trialDff = cellfun(@(x, y) (x - y) ./ y, currData.rawFl, ...
        mat2cell(currData.trialBaseline, ones(size(currData.trialBaseline))), 'uniformoutput', 0);
currData.expDff = cellfun(@(x, y) (x - y) ./ y, currData.rawFl, ...
        mat2cell(currData.expBaseline, ones(size(currData.expBaseline))), 'uniformoutput', 0);
    
wedgeData = [];
for iExp = 1:size(expMd, 1)
    currExpData = currData(strcmp(currData.expID, expMd.expID{iExp}), :);
    for iTrial = 1:max(currExpData.trialNum)
        
        % Consolidate Fl data into matrices
        td = sortrows(currExpData(currExpData.trialNum == iTrial, :), 'roiName');
        rawFlMat = cell2mat(td.rawFl');
        trialDffMat = cell2mat(td.trialDff');
        expDffMat = cell2mat(td.expDff');
        
        % Calculate bump population vector average
        angleMat = repmat((-7 * pi/8):pi/4:(7 * pi/8), size(expDffMat, 1), 1);
        [x, y] = pol2cart(angleMat, expDffMat); % Convert to cartesian coordinates to add vectors
        [theta, ~] = cart2pol(sum(x, 2), sum(y, 2));
        theta = -theta; % Inverting sign so it matches other data in plots 
        
        % Compile into table and append to wedgeData
        newRow = table(expMd.expID(iExp), iTrial, {rawFlMat}, {trialDffMat}, {expDffMat}, {theta}, ...
                {theta / pi * 4 + 4.5}, 'variableNames', {'expID', 'trialNum', 'rawFl', 'trialDff', ...
                'expDff', 'pvaRad', 'pvaWedge'});
        wedgeData = [wedgeData; newRow];
        
    end
end

% ---------- Consolidate PB glomerulus data ----------

% Calc dF/F for PB glom ROIs
currData = roiData(cellfun(@isempty, regexp(roiData.roiName, 'EB-', 'match')), :);
currData.trialDff = cellfun(@(x, y) (x - y) ./ y, currData.rawFl, ...
        mat2cell(currData.trialBaseline, ones(size(currData.trialBaseline))), 'uniformoutput', 0);
currData.expDff = cellfun(@(x, y) (x - y) ./ y, currData.rawFl, ...
        mat2cell(currData.expBaseline, ones(size(currData.expBaseline))), 'uniformoutput', 0);
    
glomData = [];
for iExp = 1:size(expMd, 1)
    currExpData = currData(strcmp(currData.expID, expMd.expID{iExp}), :);
    for iTrial = 1:max(currExpData.trialNum)
        
        % Consolidate Fl data into matrices
        td = sortrows(currExpData(currExpData.trialNum == iTrial, :), 'roiName');
        rawFlMat = cell2mat(td.rawFl');
        trialDffMat = cell2mat(td.trialDff');
        expDffMat = cell2mat(td.expDff');
        
        % Sort into paired order
        rawFlMat = rawFlMat(:, [1:9, 16:-1:10]);
        trialDffMat = trialDffMat(:, [1:9, 16:-1:10]);
        expDffMat = expDffMat(:, [1:9, 16:-1:10]);
        
        % Compile into table and append to glomData
        newRow = table(expMd.expID(iExp), iTrial, {rawFlMat}, {trialDffMat}, {expDffMat}, ...
                'variableNames', {'expID', 'trialNum', 'rawFl', 'trialDff', 'expDff'});
        glomData = [glomData; newRow];
        
    end
end
disp('Loading complete')

end

