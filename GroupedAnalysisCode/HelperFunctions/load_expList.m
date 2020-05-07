function expList = load_expList(varargin)
% Loads an expList .csv file and optionally filters by a 'groupName' column
% 
% Inputs (all optional):
%   'parentDir' (default: 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData')
%   'fileName' (default: 'fullExpList.csv')
%   'groupName' (default: no filtering)
%
p = inputParser;
addParameter(p, 'parentDir', 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData');
addParameter(p, 'fileName', 'fullExpList.csv');
addParameter(p, 'groupName', '');
parse(p, varargin{:});
parentDir = p.Results.parentDir;
fileName = p.Results.fileName;
groupName = p.Results.groupName;


expList = readtable(fullfile(parentDir, fileName), 'delimiter', ',');
if ~isempty(groupName)
   expList = expList(strcmp(expList.groupName, groupName), :); 
end



end