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
behavAnnot(end + 1, 1:11) = 1:11; % To keep colors consistent
figure(1); clf; imagesc(behavAnnot);
cMap = [rgb('Indigo'); rgb('Gold');rgb('DarkRed');rgb('Cyan');rgb('Yellow');rgb('Green');rgb('Red');rgb('blue');rgb('Gold');rgb('black')];
colormap(cMap);
yl = ylim(); ylim(yl - [0 1]) % To cut off the extra row that's just there for colormapping

% Display available annotation event types
disp(annotationTypeSummary)

%%

% Choose active alignment events
activeEventTypes = [4];
alignEventSummary = annotationTypeSummary(activeEventTypes, 2);
disp(alignEventSummary)

% Choose active filter events
activeFilterTypes = [2 3];
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
analysisWindows = []; overshoots = []; filterMatches = [];


% % Closed Loop stim
% analysisWindows(end+1,:) = [ 1  2 ];
% overshoots(end+1)        = 1;
% filterMatches{end+1} = [];
% 
% % No Stim
% analysisWindows(end+1,:) = [ 2  3 ];
% overshoots(end+1)        = 1;
% filterMatches{end+1} = [];

% % Odor A
% analysisWindows(end+1,:) = [ 2  3 ];
% overshoots(end+1)        = 1;
% filterMatches{end+1} = [];
% 
% % Odor B
% analysisWindows(end+1,:) = [ 2  3 ];
% overshoots(end+1)        = 1;
% filterMatches{end+1} = [];
% % % 
% 
% % Sound
% analysisWindows(end+1,:) = [ 2  3 ];
% overshoots(end+1)        = 1;
% filterMatches{end+1} = [];

% % All behavior
% analysisWindows(end+1,:) = [ 1  2 ];
% overshoots(end+1)        = 0;
% filterMatches{end+1} = [];
% 
% Locomotion
analysisWindows(end+1,:) = [ 2  3 ];
overshoots(end+1)        = 0;
filterMatches{end+1} = [];
% % 
% % Isolated Movement
% analysisWindows(end+1,:) = [ 2  2 ];
% overshoots(end+1)        = 0;
% filterMatches{end+1} = [];

% % Grooming
% analysisWindows(end+1,:) = [ 2  3 ];
% overshoots(end+1)        = 0;
% filterMatches{end+1} = [];

% ==================================================================================================
% ------------- FILTER EVENTS -----------------
% ==================================================================================================

% Create filters for different condition components
allFilts = []; allFiltNames = []; filtWindows = [];
% 
% % Closed Loop stim
% withStim =  [ 0  1  0 ];
% noStim =    [-1 -1  0 ];
% anyStim =   [ 0  0  0 ];
% filtWindows(end+1,:)     = [  0  0  ];
% allFilts{end+1} = [withStim; noStim]; %
% allFiltNames{end+1} = {'WithStim', 'NoStim'}; %

% % Odor
% withOdor =  [ 0  1  0 ];
% noOdor =    [-1 -1  0 ];
% anyOdor =   [ 0  0  0 ];
% filtWindows(end+1,:)     = [  0  0  ];
% allFilts{end+1} = [withOdor; noOdor]; %
% allFiltNames{end+1} = {'WithOdor', 'NoOdor'}; %
% 
% Odor A
withOdorA =  [ 0  1  0 ];
noOdorA =    [-1 -1  0 ];
anyOdorA =   [ 0  0  0 ];
filtWindows(end+1,:)     = [  0  0  ];
allFilts{end+1} = [withOdorA; noOdorA]; %
allFiltNames{end+1} = {'WithOdorA', 'NoOdorA'}; %
% 
% Odor B
withOdorB =  [ 0  1  0 ];
noOdorB =    [-1 -1  0 ];
anyOdorB =   [ 0  0  0 ];
filtWindows(end+1,:)     = [  0  0  ];
allFilts{end+1} = [withOdorB; noOdorB]; %
allFiltNames{end+1} = {'WithOdorB', 'NoOdorB'}; %
% % % % 
% % Sound
% withSound =  [ 0  1  0 ];
% noSound =    [-1 -1  0 ];
% anySound =   [ 0  0  0 ];
% filtWindows(end+1,:)     = [  0  0  ];
% allFilts{end+1} = [withSound; noSound]; %
% allFiltNames{end+1} = {'WithSound', 'NoSound'}; %

% % No stim
% withNoStim =  [ 0  1  0 ];
% noNoStim =    [-1 -1  0 ];
% anyNoStim =   [ 0  0  0 ];
% filtWindows(end+1,:)     = [  0  0  ];
% allFilts{end+1} = [withNoStim; noNoStim]; %
% allFiltNames{end+1} = {'WithNoStim', 'NoNoStim'}; %
% 
% % Locomotion
% startLoc = [-1  1  0 ];
% endLoc =   [ 1  0 -1 ];
% contLoc =  [ 1  1  0 ];
% noLoc =    [-1 -1 -1 ];
% anyLoc =   [ 0  0  0 ];
% withLoc =  [ 0  1  0 ];
% filtWindows(end+1,:)     = [ 0  0 ];
% allFilts{end+1} = [noLoc; startLoc; contLoc]; %endMove; anyMove;startLoc; 
% allFiltNames{end+1} = { 'NoLoc', 'startLoc', 'contLoc'}; %'EndMove', 'AnyMove','StartLoc',

% % Isolated movement
% startIsoMove = [-1  1  0 ];
% endIsoMove =   [ 1  0 -1 ];
% contIsoMove =  [ 1  1  0 ];
% noIsoMove =    [-1 -1 -1 ];
% anyIsoMove =   [ 0  0  0 ];
% withIsoMove =  [ 0  1  0 ];
% afterIsoMove = [ 1  0  0 ];
% beforeIsoMove =[-1  1  0 ];
% filtWindows(end+1,:)     = [  1  1  ];
% allFilts{end+1} = [noIsoMove; beforeIsoMove; afterIsoMove]; %endMove; anyMove;
% allFiltNames{end+1} = {'NoIsoMove', 'beforeIsoMove', 'afterIsoMove'}; %'EndMove', , 'AnyMove'
% 
% % Grooming
% startGroom = [-1  1  0 ];
% endGroom =   [ 1  0 -1 ];
% contGroom =  [ 1  1  0 ];
% noGroom =    [-1 -1 -1 ];
% anyGroom =   [ 0  0  0 ];
% withGroom =  [ 0  1  0 ];
% filtWindows(end+1,:)     = [  1  1  ];
% allFilts{end+1} = [noGroom; withGroom]; %endMove; anyMove;
% allFiltNames{end+1} = {'NoGroom', 'withGroom'}; %'EndMove', , 'AnyMove'
% 
% % 
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
