function allRoiData = load_roi_data(expList, parentDir)

if nargin < 2
    parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\all_experiments';
end
if isa(expList, 'table')
    expList = expList.expID;
end

allRoiData = [];
for iExp = 1:size(expList, 1)
   currExpID = expList{iExp};
   roiDataFile = fullfile(parentDir, [currExpID, '_roiData.mat']);
   if exist(roiDataFile, 'file')
       disp(['Loading ', currExpID, '...'])
       load(roiDataFile, 'roiData');
       allRoiData = [allRoiData; roiData];
   else
       disp(['Skipping ', currExpID, '...file not found']);
   end
end

end