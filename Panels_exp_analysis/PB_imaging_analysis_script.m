%% LOAD DATA
% Load all analysis data for one or more experiments, and generate arrays of fluorescence data for 
% ordered EB wedges and individual PB glomeruli. 

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData_60D05_7f';
expList = {'20210118-1', '20210118-2'};
figDir = fullfile(parentDir, 'Figs');

[expMd, trialMd, roiData, ftData, flailingEvents, panelsMetadata, wedgeData, glomData] = ...
        load_PB_data(parentDir, expList);

% Load file with info about the details of the bath applications
drugTimingMd = readtable(fullfile(parentDir, 'drug_timing_data.csv'), 'delimiter', ',');

%---------------------------------------------------------------------------------------------------

%% CHECK CORRELATION BETWEEN GLOMERULI
% Make sure the glomeruli were identified correctly by looking at the correlation between the
% left and right matching glomeruli.
%---------------------------------------------------------------------------------------------------
expList = unique(expMd.expID);%{'20201210-1'}%

for iExp = 1:numel(expList)
    check_glom_correlation(expList{iExp}, expMd, roiData);
end 

%% PLOT OVERVIEW OF VISUAL TUNING AND MOVEMENT FOR A SINGLE TRIAL


expID = '20210118-1';
trialNum = 8;
sourceData = wedgeData; % glomData or wedgeData
p = [];
p.flType = 'expDff';  % rawFl, trialDff, or expDff
p.plotPVA = 1;
p.plotMeanPVA = 1;
p.useFlow = 1;
p.flMax = [];
p.smWin = 5;
p.figNum = [];

%---------------------------------------------------------------------------------------------------

currData = inner_join(sourceData, expMd, trialMd, panelsMetadata);
currData = outerjoin(currData, ftData, 'type', 'left', 'mergekeys', 1);
trialData = currData(strcmp(currData.expID, expID) & currExpData.trialNum == trialNum, :);
plot_single_trial_visual_tuning_summary(trialData, p);








