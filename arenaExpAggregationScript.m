%% Get mean walking speed and locomotion liklihood for a whole list of experiments

parentDir = 'D:\Dropbox (HMS)\Behavior_Experiments';

locThresh = 1;
smWin = 9;

% expDates = ["2019_03_07_exp_4", ...
%     "2019_03_11_exp_4", ...
%     "2019_03_22_exp_1", ...
%     "2019_03_24_exp_1", ...
%     "2019_03_26_exp_3", ...
%     ];

expDates = ["2019_03_14_exp_1", ...
        "2019_03_15_exp_1", ...
        "2019_03_15_exp_2", ...
        "2019_03_20_exp_1", ...
        "2019_03_20_exp_2", ...
        "2019_03_21_exp_1", ...
        "2019_03_21_exp_2", ...
        "2019_03_21_exp_2", ...
        "2019_03_24_exp_3", ...
        "2019_03_26_exp_4", ...
        "2019_03_27_exp_1", ...
        "2019_03_27_exp_3", ...
    ];

% allTids = {[0 1], ...
%     [0 1], ...
%     [0 1], ...
%     [0 1], ...
%     [0 1], ...
%     };

allTids = {[0], ...
    [0 1 2], ...
    [0 1 2], ...
    [0 1 2], ...
    [0 1], ...
    [1 2], ...
    [1], ...
    0, ...
    0, ...
    [0 1], ...
    1, ...
    1, ...
    };

% try
avgSpeedBaseline = []; avgSpeedStim = []; avgLocBaseline = []; avgLocStim = [];
for iExp = 1:numel(expDates)
    expDate = expDates{iExp};
    tids = allTids{iExp};
    trxData = [];
    for iTrial = 1:numel(tids)
        tid = tids(iTrial);
        
        % Load trial metadata
        metadata = load(fullfile(parentDir, expDate, ['tid_', num2str(tid), '_metadata.mat']));
        FRAME_RATE = metadata.FRAME_RATE;
        stimOnsets = metadata.stimOnsets;
        stimDurs = metadata.stimDurs;
        if ~isempty(stimDurs)
            stimDur = stimDurs(1);
        end
        
        % Load CTrax output file
        if exist(fullfile(parentDir, expDate, ['tid_', num2str(tid), '_with_stats.mat']), 'file')
            [trx, ~ , ~] = load_tracks(fullfile(parentDir, expDate, ['tid_', num2str(tid), '_with_stats.mat'])); % --> trx
        elseif exist(fullfile(parentDir, expDate, ['fixed_tid_', num2str(tid), '.mat']), 'file')
            [trx, ~ , ~] = load_tracks(fullfile(parentDir, expDate, ['fixed_tid_', num2str(tid), '.mat'])); % --> trx
        else
            [trx, ~ , ~] = load_tracks(fullfile(parentDir, expDate, ['tid_', num2str(tid), '.mat'])); % --> trx
        end
        
        trxData = [trxData, trx];
        disp(['nFlies = ', num2str(numel(trx))]);
    end
    
    % EXTRACT DATA
    flyBirths = [trxData.firstframe];
    flyDeaths = [trxData.endframe];
    nFrames = max(flyDeaths);
    stimOnsetFrames = round(stimOnsets * FRAME_RATE);
    stimOffsetFrames = round((stimOnsets + stimDur) * FRAME_RATE);
    stimOffsetFrames(stimOffsetFrames > nFrames) = nFrames;
    nStims = numel(stimOnsets);
    stimOnFrames = zeros(nFrames, 1);
    for iStim = 1:nStims
        stimOnFrames(stimOnsetFrames(iStim):stimOffsetFrames(iStim) - 1) = 1;
    end
    
    % Figure out which flies are tracked throughout entire trial
    goodFlies = flyBirths == 1 & flyDeaths >= (nFrames - 10); % Allow up to 10 dropped frames
    goodFlies = logical(goodFlies);
    goodFlyData = trxData(goodFlies);
    nGoodFrames = min(flyDeaths(goodFlies));
    disp('-------------------------------------------')
    disp(expDate);
    disp(['Number of good fly tracks: ', num2str(sum(goodFlies))]);
    disp(['Number of bad fly tracks: ', num2str(sum(~goodFlies))]);
    
    % Extract speed and locomotion data
    speedData = []; locData = [];
    for iFly = 1:sum(goodFlies)
        speedData(:, iFly) = goodFlyData(iFly).velmag_ctr(1:nGoodFrames-1); % --> [frame, fly]
    end
    speedDataBaseline = speedData;
    speedDataStim = speedData;
    speedDataBaseline(logical(stimOnFrames(1:size(speedData, 1))), :) = nan;
    speedDataStim(~logical(stimOnFrames(1:size(speedData, 1))), :) = nan;
    
    smSpeedBaseline = movmean(speedDataBaseline, smWin, 1, 'omitnan');
    smSpeedStim = movmean(speedDataStim, smWin, 1, 'omitnan');
    locDataBaseline = smSpeedBaseline > locThresh;
    locDataStim = smSpeedStim > locThresh;
    
    avgSpeedBaseline{iExp} = mean(smSpeedBaseline, 1, 'omitnan');
    avgSpeedStim{iExp} = mean(smSpeedStim, 1, 'omitnan');
    avgLocBaseline{iExp} = mean(locDataBaseline, 1, 'omitnan');
    avgLocStim{iExp} = mean(locDataStim, 1, 'omitnan');
end
% catch foldME; rethrow(foldME); end

%%
linMeanSpeed = []; linMeanLoc = []; expGroupVars = []; varNameVars = [];
for iExp = 1:numel(expDates)
   currSpeedData = [avgSpeedBaseline{iExp}, avgSpeedStim{iExp}]';
   
   expGroup = repmat({expDates{iExp}}, 1, numel(currSpeedData));
   varNameGroups = [repmat({'Baseline'}, 1, numel(avgSpeedBaseline{iExp})), ...
                    repmat({'Stim'}, 1, numel(avgSpeedStim{iExp}))];
                
   
   linMeanSpeed = [linMeanSpeed; currSpeedData];
   expGroupVars = [expGroupVars, expGroup];
   varNameVars = [varNameVars, varNameGroups];   
end

% groups = {regexprep(expGroupVars, '(\_|2019\_0)', '')';varNameVars'};



%%
speedData = []; catIdx = {};
for iExp = 1:numel(expDates)
   speedData{iExp} = [avgLocBaseline{iExp}, avgLocStim{iExp}];
   catIdx = [catIdx, repmat({'Baseline frames'}, 1, numel(avgSpeedBaseline{iExp})), ...
                      repmat({'Stim frames'}, 1, numel(avgSpeedStim{iExp}))]; 
end

distIdx = [ones(1,32), 2*ones(1,20), 3*ones(1,20), 4*ones(1,28), 5*ones(1,28)];
 

f = figure(1);clf;
f.Color = [1 1 1];
handles = plotSpread(speedData, 'categoryIdx', catIdx, 'xNames', ...
        regexprep(expDates, {'2019\_', '\_'}, {'', '\.'}), 'categoryMarkers', {'o', '*'}, ...
        'categoryColors', {'b', rgb('green')}, 'showMM', 4);
legend('Baseline frames', 'Stim frames')
ax = handles{3};
ax.FontSize = 14;
ax.XTickLabelRotation = 90;
ax.Position(2) = ax.Position(2) * 1.2;
title('Mean locomotion speed  (PPM2 > op6s control expts)');
ylabel('Move speed (mm/sec)')

%%


speedData = []; catIdx = {};
for iExp = 1:numel(expDates)
   speedData{iExp} = [avgSpeedBaseline{iExp}];
   catIdx = [catIdx, repmat({'Baseline frames'}, 1, numel(avgSpeedBaseline{iExp})), ...
                      repmat({'Stim frames'}, 1, numel(avgSpeedStim{iExp}))]; 
end

distIdx = [ones(1,32), 2*ones(1,20), 3*ones(1,20), 4*ones(1,28), 5*ones(1,28)];
 

f = figure(1);clf;
f.Color = [1 1 1];
handles = plotSpread(speedData, 'xNames', regexprep(expDates, {'2019\_', '\_'}, {'', '\.'}), ...
         'distributionMarkers', 'o', 'showMM', 4);
% legend('Baseline frames', 'Stim frames')
ax = handles{3};
ax.FontSize = 14;
ax.XTickLabelRotation = 90;
ax.Position(2) = ax.Position(2) * 1.2;
title('Mean movement speed  (PPM2 > csChrimson expts)');
ylabel('Movement speed (mm/sec)')

save_figure(f, parentDir, 'csChrimson_speed');






















