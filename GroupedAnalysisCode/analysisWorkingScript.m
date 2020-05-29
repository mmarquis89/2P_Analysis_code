
%% Load a base aligned data object

parentDir = 'D:\Dropbox (HMS)\2P Data\Imaging Data\GroupedAnalysisData\Saved_AlignEvent_objects';
alignEventDateStr = '20200528';

% ballstop, flailing, grooming, isolatedmovement, locomotion, odor, optostim, panelsflash, soundstim
alignEventName = 'flailing';

load(fullfile(parentDir, [alignEventDateStr, '_AlignEventObj_', alignEventName, '.mat']));

%% Generate and adjust filter structure
filterDefs = alignObj.create_filterDefs();

% General fields
filterDefs.expID = [];
filterDefs.trialNum = [];
filterDefs.expName = 'D-ANT';
filterDefs.moveSpeed = [];
filterDefs.yawSpeed = [];
filterDefs.roiName = 'ANT';

% Event filter vectors
filterDefs.ballstop = 0;
filterDefs.grooming = 0;
filterDefs.isolatedmovement = 0;
filterDefs.locomotion = 0;
filterDefs.odor = 0;
filterDefs.optostim = 0;
filterDefs.panelsflash = 0;
filterDefs.soundstim = [];

% OneNote experiment metadata
filterDefs.genotype = [];
filterDefs.ageHrs = [];
filterDefs.foodType = [];
filterDefs.bedtime = [];
filterDefs.starvationTimeHrs = [];
filterDefs.olfactometerVersion = [];
filterDefs.bathTemp = [];

% Alignment event-specific fields
if strcmp(alignObj.alignEventName, 'odor')
    filterDefs.odor = [];
    filterDefs.concentration = [];
end

%% Generate filtered and aligned event data table

analysisWin = [2 2];

dataTable = alignObj.output_analysis_subset(analysisWin, filterDefs);

%% 






