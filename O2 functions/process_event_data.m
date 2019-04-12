%% ===================================================================================================
%%% PROCESS EVENT DATA BASED ON TRIAL INTERACTIONS
%% ===================================================================================================

% Required inputs:
%       A .mat file with variables 'annotationTypes' and 'annotationTypeSummary'
%
% Outputs:
%       A .mat file with the following variables:
%           alignEventSummary
%           filterEventSummary
%           primaryEventNames
%           eventLists
%           nEventTypes
%           condNames
%           onsetFilterVecs
%           offsetFilterVecs
%===================================================================================================

% Load annotation type data
[annotTypeFileName, parentDir] = uigetfile('D:\Dropbox (HMS)\2P Data\Imaging Data', 'Select an annotationType data file');
load(fullfile(parentDir, annotTypeFileName)) % Vars 'annotationTypes', 'annotationTypeSummary'

% Load analysis metadata
if exist(fullfile(parentDir, 'analysisMetadata.mat'), 'file')
    load(fullfile(parentDir, 'analysisMetadata.mat')) % 'analysisMetadata'
else
    [analysisMetadataFile, parentDir] = uigetfile('D:\Dropbox (HMS)\2P Data\ImagingData');
    load(fullfile(parentDir, analysisMetadataFile)); % 'analysisMetadata'
end
volumeRate = analysisMetadata.volumeRate;


%%
% Show preview of behavior annotation data
behavAnnot = annotationTypes{end}.frameAnnotArr;
behavAnnot(end + 1, 1:10) = 1:10; % To keep colors consistent
figure(1); clf; imagesc(behavAnnot);
cMap = [rgb('Indigo'); rgb('orange');rgb('DarkRed');rgb('Cyan');rgb('Yellow');rgb('Green');rgb('Red');rgb('blue');rgb('Gold');rgb('black')];
colormap(cMap);
yl = ylim(); ylim(yl - [0 1]) % To cut off the extra row that's just there for colormapping

% Display available annotation event types
disp(annotationTypeSummary)

%%

% Choose active alignment events
activeEventTypes = [2];
alignEventSummary = annotationTypeSummary(activeEventTypes, 2);
disp(alignEventSummary)

% Choose active filter events
activeFilterTypes = [4];
filterEventSummary = annotationTypeSummary(activeFilterTypes, 2);
disp(filterEventSummary)

% Generate output file name 
saveFileName = 'EventData_Align';
for iAlign = 1:numel(alignEventSummary) 
    saveFileName = [saveFileName, '_', cell2mat(alignEventSummary{iAlign,1})];
end
saveFileName = [saveFileName, '_Filter'];
for iFilt = 1:numel(filterEventSummary) 
    saveFileName = [saveFileName, '_', cell2mat(filterEventSummary{iFilt,1})];
end
saveFileName = [saveFileName, '.mat'];

% ==================================================================================================
% --------- ALIGNMENT EVENTS -------------
% ==================================================================================================
analysisWindows = []; overshoots = []; filterMatches = []; alignEventInfo = [];

% Closed loop stim
alignEventInfo.ClosedLoop.analysisWindows = [1 2];
alignEventInfo.ClosedLoop.overshoots = 1;
alignEventInfo.ClosedLoop.filterMatches = {[]};

% No stim
alignEventInfo.NoStim.analysisWindows = [2 3];
alignEventInfo.NoStim.overshoots = 1;
alignEventInfo.NoStim.filterMatches = {[]};

% Odor A
alignEventInfo.OdorA.analysisWindows = [2 3];
alignEventInfo.OdorA.overshoots = 1;
alignEventInfo.OdorA.filterMatches = {[]};

% Odor AB pair
alignEventInfo.OdorABPair.analysisWindows = [2 4];
alignEventInfo.OdorABPair.overshoots = 1;
alignEventInfo.OdorABPair.filterMatches = {[]};

% Odor B
alignEventInfo.OdorB.analysisWindows = [2 3];
alignEventInfo.OdorB.overshoots = 1;
alignEventInfo.OdorB.filterMatches = {[]};

% Sound
alignEventInfo.Sound.analysisWindows = [2 3];
alignEventInfo.Sound.overshoots = 1;
alignEventInfo.Sound.filterMatches = {[]};

% All behavior
alignEventInfo.move.analysisWindows = [1 2];
alignEventInfo.move.overshoots = 0;
alignEventInfo.move.filterMatches = {[]};

% Locomotion
alignEventInfo.locomotion.analysisWindows = [2 3];
alignEventInfo.locomotion.overshoots = 0;
alignEventInfo.locomotion.filterMatches = {[]};

% Isolated movement
alignEventInfo.isoMove.analysisWindows = [2 3];
alignEventInfo.isoMove.overshoots = 0;
alignEventInfo.isoMove.filterMatches = {[]};


% Compile active event info
for iEvent = 1:numel(activeEventTypes)
    currEventName = annotationTypeSummary.AnnotationType{activeEventTypes(iEvent)};
    analysisWindows(end + 1, :) = alignEventInfo.(currEventName).analysisWindows;
    overshoots(end + 1) = alignEventInfo.(currEventName).overshoots;
    filterMatches{end + 1} = alignEventInfo.(currEventName).filterMatches{:};
end

% ==================================================================================================
% ------------- FILTER EVENTS -----------------
% ==================================================================================================

% Create filters for different condition components
allFilts = []; allFiltNames = []; filtWindows = []; filterEventInfo = [];

% Closed Loop stim
withStim =  [ 0  1  0 ];
noStim =    [-1 -1  0 ];
anyStim =   [ 0  0  0 ];
filterEventInfo.ClosedLoop.filtWindows = [  0  0  ];
filterEventInfo.ClosedLoop.allFilts = [withStim; noStim];
filterEventInfo.ClosedLoop.allFiltNames = {'WithStim', 'NoStim'};

% Odor A
withOdorA =  [ 0  1  0 ];
noOdorA =    [-1 -1  0 ];
anyOdorA =   [ 0  0  0 ];
filterEventInfo.OdorA.filtWindows = [  0  0  ];
filterEventInfo.OdorA.allFilts = [withOdorA; noOdorA]; 
filterEventInfo.OdorA.allFiltNames = {'WithOdorA', 'NoOdorA'}; 

% Odor AB Pair
withOdorABPair =  [ 0  1  0 ];
noOdorABPair =    [-1 -1  0 ];
anyOdorABPair =   [ 0  0  0 ];
filterEventInfo.OdorABPair.filtWindows = [  0  0  ];
filterEventInfo.OdorABPair.allFilts = [withOdorABPair; noOdorABPair]; 
filterEventInfo.OdorABPair.allFiltNames = {'WithOdorABPair', 'NoOdorABPair'}; 

% Odor B
withOdorB =  [ 0  1  0 ];
noOdorB =    [-1 -1  0 ];
anyOdorB =   [ 0  0  0 ];
filterEventInfo.OdorB.filtWindows = [  0  0  ];
filterEventInfo.OdorB.allFilts = [withOdorB; noOdorB]; 
filterEventInfo.OdorB.allFiltNames = {'WithOdorB', 'NoOdorB'}; 

% Sound
withSound =  [ 0  1  0 ];
noSound =    [-1 -1  0 ];
anySound =   [ 0  0  0 ];
filterEventInfo.Sound.filtWindows = [  0  0  ];
filterEventInfo.Sound.allFilts = [withSound; noSound]; 
filterEventInfo.Sound.allFiltNames = {'WithSound', 'NoSound'}; 

% No stim
withNoStim =  [ 0  1  0 ];
noNoStim =    [-1 -1  0 ];
anyNoStim =   [ 0  0  0 ];
filterEventInfo.NoStim.filtWindows = [  0  0  ];
filterEventInfo.NoStim.allFilts = [withNoStim; noNoStim]; 
filterEventInfo.NoStim.allFiltNames = {'WithNoStim', 'NoNoStim'}; 

% Locomotion
startLoc = [-1  1  0 ];
endLoc =   [ 1  0 -1 ];
contLoc =  [ 1  1  0 ];
noLoc =    [-1 -1 -1 ];
anyLoc =   [ 0  0  0 ];
withLoc =  [ 0  1  0 ];
filterEventInfo.locomotion.filtWindows = [ 1  1 ];
filterEventInfo.locomotion.allFilts = [noLoc; startLoc; contLoc]; %endMove; anyMove;startLoc; 
filterEventInfo.locomotion.allFiltNames = { 'NoLoc', 'startLoc', 'contLoc'}; %'EndMove', 'AnyMove','StartLoc',

% Isolated movement
startIsoMove = [-1  1  0 ];
endIsoMove =   [ 1  0 -1 ];
contIsoMove =  [ 1  1  0 ];
noIsoMove =    [-1 -1 -1 ];
anyIsoMove =   [ 0  0  0 ];
withIsoMove =  [ 0  1  0 ];
afterIsoMove = [ 1  0  0 ];
beforeIsoMove =[-1  1  0 ];
filterEventInfo.isoMove.filtWindows = [ 1  1 ];
filterEventInfo.isoMove.allFilts = [noIsoMove; beforeIsoMove; afterIsoMove]; %endMove; anyMove;
filterEventInfo.isoMove.allFiltNames = {'NoIsoMove', 'beforeIsoMove', 'afterIsoMove'}; %'EndMove', , 'AnyMove'

% All behavior
startMove = [-1  1  0 ];
endMove =   [ 1  0 -1 ];
contMove =  [ 1  1  0 ];
noMove =    [-1 -1 -1 ];
anyMove =   [ 0  0  0 ];
filterEventInfo.move.filtWindows = [ 1  1 ];
filterEventInfo.move.allFilts = [noMove; startMove; contMove]; %endMove; anyMove;
filterEventInfo.move.allFiltNames = {'NoMove', 'StartMove', 'ContMove'}; %'EndMove', , 'AnyMove'


% Compile active filter info
for iEvent = 1:numel(activeFilterTypes)
    currEventName = annotationTypeSummary.AnnotationType{activeFilterTypes(iEvent)};
    filtWindows(end + 1, :) = filterEventInfo.(currEventName).filtWindows;
    allFilts{end + 1} = filterEventInfo.(currEventName).allFilts;
    allFiltNames{end + 1} = filterEventInfo.(currEventName).allFiltNames;
end



% % All behavior
% startMove = [-1  1  0 ];
% endMove =   [ 1  0 -1 ];
% contMove =  [ 1  1  0 ];
% noMove =    [-1 -1 -1 ];
% anyMove =   [ 0  0  0 ];
% filtWindows(end+1,:)     = [ 0  0 ];
% allFilts{end+1} = [noMove; startMove; contMove]; %endMove; anyMove;
% allFiltNames{end+1} = {'NoMove', 'StartMove', 'ContMove'}; %'EndMove', , 'AnyMove'

%===================================================================================================
%
%                                                 |------------event------------|          
%       <--------fW(1)--------><----aW(1)---->[alignVol]<----aW(2)----><----fW(2)---->
%       [------------------fD(1)-----------------][-------fD(2)-------][----fD(3)----]
%
%===================================================================================================
    
% Get event vols for each filter type
filterEventVols = [];
for iType = 1:numel(annotationTypes)
    filterEventVols(:,:,iType) = annotationTypes{iType}.eventVols;    
end

% Compile annotation info for active events
primaryEventNames = [];
for iType = 1:numel(activeEventTypes)
    currType = activeEventTypes(iType);
    eventLists{iType} = annotationTypes{currType}.eventList;
    primaryEventNames{iType} = annotationTypes{currType}.name;
end

nEventTypes = numel(primaryEventNames);
allCondFilters = []; condNames = []; activeFilterEventVols = []; onsetFilterVecs = []; offsetFilterVecs = [];
for iType = 1:nEventTypes
    
    % Eliminate filter types that overlap with the current primary event
    matchedFilterTypes = zeros(1, numel(allFilts));
    matchedFilterTypes(filterMatches{iType}) = 1;
    currFilts = allFilts(~matchedFilterTypes);
    currFiltNames = allFiltNames(~matchedFilterTypes);
    
    % Combine primary events and filter types into all filter conditions
    [allCondFilters{iType}, condNames{iType}] =  create_filter_conditions(primaryEventNames{iType}, currFilts, currFiltNames);
    
    % Select the correct eventVols for each primary event type
    activeFilterEventVols = filterEventVols(:,:,activeFilterTypes(~matchedFilterTypes));    
    
    % Get onset filter vecs for each condition
    for iCond = 1:numel(condNames{iType})
        onsetFilterVecs{iType}(:, iCond) = filter_event_data(eventLists{iType}, activeFilterEventVols, analysisWindows(iType, :), ...
            filtWindows(~matchedFilterTypes, :), allCondFilters{iType}{iCond}, 'volumeRate', volumeRate, 'overshoot', overshoots(iType));
    end
    
    % Get offset filter vecs for each condition
    for iCond = 1:numel(condNames{iType})
        offsetFilterVecs{iType}(:, iCond) = filter_event_data(eventLists{iType}, activeFilterEventVols, analysisWindows(iType, :), ...
            filtWindows(~matchedFilterTypes, :), allCondFilters{iType}{iCond}, 'volumeRate', volumeRate, 'overshoot', overshoots(iType), 'offsetAlign', 1);
    end
end

% Save data and metadata in a .mat file
save(fullfile(parentDir, saveFileName), 'alignEventSummary', 'filterEventSummary', 'primaryEventNames', 'eventLists', ...
    'nEventTypes', 'condNames', 'onsetFilterVecs', 'offsetFilterVecs', 'analysisWindows', 'filtWindows', '-v7.3');

clear withOdor noOdor anyOdor startMove endMove contMove noMove anyMove withLaser noLaser anyLaser currType matchedFilterTypes currFilts currFiltNames



