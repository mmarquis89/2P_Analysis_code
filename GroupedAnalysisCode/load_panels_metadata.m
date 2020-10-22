function panelsMetadata = load_panels_metadata(expList, parentDir)
% ==================================================================================================   
%  Loads panels metadata files for a set of expIDs and compiles them into a single table. 
%
%  INPUTS: 
%       expList                 = cell array of expIDs (YYYYMMDD-expNum format) to load data for.
%                                 Can provide a table with a field named "expID" instead.
%
%       parentDir (optional)    = override default location for source files
%
%
%  OUTPUTS:
%       panelsMetadata = table with these columns of FicTrac data for each trial:
%                                   expID
%                                   trialNum
%                                   panelsMode
%                                   pattern
%                                   xDimPosFunc
%                                   yDimPosFunc
%                                   panelsPosX
%                                   panelsPosY
%
% ==================================================================================================

% Parse/process input arguments
if nargin < 2
    parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';
end
if isa(expList, 'table')
    expList = expList.expID;
end

% Load panels metadata file for each expID if it exists
panelsMd = [];
disp('------------------------------------------');
disp('Loading panels metadata files...')
for iExp = 1:numel(expList)
   currExpID = expList{iExp};
   panelsMdFile = fullfile(parentDir, [currExpID, '_panelsMetadata.mat']);
   if exist(panelsMdFile, 'file')
       disp(['Loading ', currExpID, '...'])
       load(panelsMdFile, 'panelsMetadata');
       panelsMd = [panelsMd; panelsMetadata];
   else
       disp(['Skipping ', currExpID, '...file not found']);
   end
end
disp('Panels metadata loaded')
end