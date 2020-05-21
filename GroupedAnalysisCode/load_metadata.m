function [expMetadata, trialMetadata] =  load_metadata(expList, parentDir)
% ==================================================================================================   
%  Loads the experiment and trial metadata files for a set of expIDs and compiles them into two 
%  tables.
%
%  INPUTS: 
%       expList                 = a cell array of expIDs (YYYYMMDD-expNum format) to load data for.
%                                 Can provide a table with a field named "expID" instead.
%
%       parentDir (optional)    = override default location for source files
%
%
%  OUTPUTS:
%       expMetadata             = table with these columns of experiment-specific info:
%                                       expID
%                                       expName
%                                       daqSampRate
%                                       panelsDisplayRate
%                                       volumeRate
%                                       nPlanes
%                                       nTrials
%
%       trialMetadata           = table with these columns of trial-specific info:
%                                       expID
%                                       trialNum
%                                       trialDuration
%                                       nVolumes
%                                       nDaqSamples
%                                       nPanelsFrames
%                                       usingOptoStim
%                                       usingPanels
%                                       using2P
%                                       originalTrialCount
%                                       pmtShutoffVols
%
% ==================================================================================================

if nargin < 2
   parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments'; 
end
if isa(expList, 'table')
   expList = expList.expID; 
end

expMetadata = [];
trialMetadata = [];
for iExp = 1:size(expList, 1)
   currExpID = expList{iExp};
   expMdFile = fullfile(parentDir, [currExpID, '_expMetadata.csv']);
   if exist(expMdFile, 'file')
       disp(['Adding ', currExpID, '...'])
       currExpMd = readtable(expMdFile, 'delimiter', ',');
       expMetadata = [expMetadata; currExpMd];
       
       trialMdFile = fullfile(parentDir, [currExpID, '_trialMetadata.csv']);
       currTrialMd = readtable(trialMdFile, 'delimiter', ',');
       trialMetadata = [trialMetadata; currTrialMd];
   else
       disp(['Skipping ', currExpID, '...file not found']);
   end
end

end