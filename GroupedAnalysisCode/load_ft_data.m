function allFtData = load_ft_data(expList, parentDir)
% ==================================================================================================   
%  Loads FicTrac data files for a set of expIDs and compiles them into a single table. 
%
%  INPUTS: 
%       expList                 = cell array of expIDs (YYYYMMDD-expNum format) to load data for.
%                                 Can provide a table with a field named "expID" instead.
%
%       parentDir (optional)    = override default location for source files
%
%
%  OUTPUTS:
%       allFtData               = table with these columns of FicTrac data for each trial:
%                                   expID
%                                   trialNum
%                                   intX
%                                   intY
%                                   intHD
%                                   moveSpeed
%                                   intFwMove
%                                   intSideMove
%                                   yawSpeed
%                                   fwSpeed
%                                   sideSpeed
%                                   frameTimes   
%                                   badVidFrames    
%                                   meanFlow                                       
%
% ==================================================================================================

% Parse/process input arguments
if nargin < 2
    parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';
end
if isa(expList, 'table')
    expList = expList.expID;
end

% Load FicTrac data file for each expID if it exists
allFtData = [];
disp('------------------------------------------');
disp('Loading FicTrac data files...')
for iExp = 1:numel(expList)
   currExpID = expList{iExp};
   ftDataFile = fullfile(parentDir, [currExpID, '_ficTracData.mat']);
   if exist(ftDataFile, 'file')
       disp(['Loading ', currExpID, '...'])
       load(ftDataFile, 'ftData');
       allFtData = [allFtData; ftData];
   else
       disp(['Skipping ', currExpID, '...file not found']);
   end
end
disp('FicTac data loaded')
end