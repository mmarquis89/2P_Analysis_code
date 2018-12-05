%%% Closed loop analyis testing code
% 
% fictracData --> [???, odorA valve, odorB shutoff, NO valve, yaw]
% Arduino sampling rate = 4000 Hz;
%
% OdorOn voltage ranges: 
%       2018_07_07_exp_1_tid_1-2:   [1.77 3.6] V 
%       2018_07_07_exp_1_tid_3:     [1.77 5.2] V
%       2018_07_07_exp_1_tid_4:     [0 3.5] V
%       2018_07_07_exp_2:           [0 5] V

odorOnRange = [0 5];
tid = 1;

% Load Fictrac data from behavior vids
ftVidData = load_fictrac_data();

% Load data
% parentDir = 'B:\Dropbox (HMS)\2P Data\Behavior Vids\2018_07_07_exp_1';
% dataFile = 'fictracData_20180707_143248_sid_0_tid_1_Closed-Loop-Odor-A.mat';
% load(fullfile(parentDir, dataFile)); % 'fictracData'
ftData = fictracData;

%%

% Set constants
FT_SAMP_RATE = 4000;
FRAME_RATE = 25;

% Downsample closed loop data to match behavior vids
dsFactor = round(FT_SAMP_RATE / FRAME_RATE);
ftDataDS = ftData(1:dsFactor:end, :);

% Initial processing steps
odorOnRangeRad = odorOnRange * (1/max(ftData(:,5))) * 2*pi;
radYaw = ftDataDS(:,5) * (1/max(ftDataDS(:,5))) * 2*pi;
odorOn = round(ftDataDS(:,4))';

% smYaw = smooth(unwrap(radYaw/max(radYaw))*2, ftDataSampRate/10)';

% Get onset/offset/event times
odorEvents = logical(odorOn);
odorEvents(1) = 0;
odorEvents(end) = 0;
odorOnStr = num2str(odorEvents);
odorOnStr = odorOnStr(~isspace(odorOnStr));
[onsetInds, offsetInds] = regexp(odorOnStr, '01+0');
odorOnsets = zeros(size(odorOn));
odorOffsets = zeros(size(odorOn));
odorOnsets(onsetInds + 1) = 1;
odorOffsets(offsetInds) = 1;

% Make event list of odor presentation
odorEventList = create_event_list(odorOnsets, odorOffsets);
eventDurs = odorEventList(:,2) - odorEventList(:,1);
eventDursSec = eventDurs / FRAME_RATE;

% Plot yaw
f = figure(1); clf; 
subaxis(2,1,1); hold on; 
plot(radYaw, '-.');
ax = gca();
ax.XTick = 1:FRAME_RATE:numel(odorOn);
ax.XTickLabel = (1:1:ceil((numel(odorOn) / FRAME_RATE)));

% Shade around odor stims
uwYaw = unwrap(radYaw);
for iEvent = 1:size(odorEventList, 1)
    startSamp = odorEventList(iEvent, 1);
    endSamp = odorEventList(iEvent, 2);
    xData = [startSamp, startSamp, endSamp, endSamp];
    yData = [odorOnRangeRad(1), odorOnRangeRad(2), odorOnRangeRad(2), odorOnRangeRad(1)];
    fill(xData, yData, rgb('red'), 'facealpha', 0.5, 'edgealpha', 0)
end

subaxis(2,1,2); hold on;
fwSmooth = smooth(ftVidData.fwSpeed(:,tid), 3);
fwNorm = fwSmooth / max(fwSmooth);
plot(fwNorm)

ax = gca();
ax.XTick = 1:FRAME_RATE:numel(odorOn);
ax.XTickLabel = (1:1:ceil((numel(odorOn) / FRAME_RATE)));

% Shade around odor stims
uwYaw = unwrap(radYaw);
for iEvent = 1:size(odorEventList, 1)
    startSamp = odorEventList(iEvent, 1);
    endSamp = odorEventList(iEvent, 2);
    xData = [startSamp, startSamp, endSamp, endSamp];
    yData = [0 1 1 0];
    fill(xData, yData, rgb('red'), 'facealpha', 0.5, 'edgealpha', 0)
end

% % Plot yaw around every event
% for iEvent = 1:size(odorEventList, 1)
% %    figure(iEvent); clf; hold on;
%    
%    
%     
% end


