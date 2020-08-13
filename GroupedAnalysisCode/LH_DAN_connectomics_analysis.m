
% CONCATENATE AND RE-SAVE CONNECTIVTY DATA FOR ALL LH-DANS

try
% Load all connectivity data
parentDir = 'C:\Users\Wilson Lab\Google Drive\Lab Work\LH-DAN_connectomics_analysis';
neuron_names = {'TypeB', 'TypeD', 'TypeF_1', 'TypeF_2'};
neuron_IDs = [328533761, 294436967, 792040520, 950229431];
tbl = [];
for iName = 1:numel(neuron_names)
    
    % Identify names of current .csv files
    curr_name = neuron_names{iName};
    us_file_name = [curr_name, '_upstream_partners.csv'];
    ds_file_name = [curr_name, '_downstream_partners.csv'];
    
    % Load upstream partners
    us_partners = readtable(fullfile(parentDir, us_file_name), 'delimiter', ',');
    ID_table = repmat(table({curr_name}, neuron_IDs(iName), {'upstream'}, 'variableNames', ...
            {'DAN_name', 'DAN_ID', 'partnerDirection'}), size(us_partners, 1), 1);
    us_partners = [ID_table, us_partners];
        
    % Load downstream partners
    ds_partners = readtable(fullfile(parentDir, ds_file_name), 'delimiter', ',');
    ID_table = repmat(table({curr_name}, neuron_IDs(iName), {'downstream'}, 'variableNames', ...
            {'DAN_name', 'DAN_ID', 'partnerDirection'}), size(ds_partners, 1), 1);
    ds_partners = [ID_table, ds_partners];
    
    % Add both to main table
    tbl = [tbl; us_partners];
    tbl = [tbl; ds_partners];
end
    
save_file = fullfile(parentDir, 'full_connectivity_table.csv');
% writetable(tbl, save_file);
catch ME; rethrow(ME);end

%% SAME INITIAL PROCESSING, BUT FOR COMPARISON CELLS

try
    
% Load all connectivity data
parentDir = 'C:\Users\Wilson Lab\Google Drive\Lab Work\LH-DAN_connectomics_analysis';
neuron_names = {'DM3', 'DM4', 'DM6', 'LHPV5a1', 'LHPV5g1', 'LHAD1b1', 'PFL1'};
neuron_IDs = {'755518957', '573333835', 'adPNs', '609941228', '728625156', '543412198', '481121605'};
% fileNamePrefixes = {'DM3_755518957', 'DM4_573333835', 'DM6_adPN'};
tbl = [];
for iName = 1:numel(neuron_names)
    
    % Identify names of current .csv files
    curr_name = neuron_names{iName};
    us_file_name = [neuron_names{iName}, '_', neuron_IDs{iName}, '_upstream_partners.csv'];
    ds_file_name = [neuron_names{iName}, '_', neuron_IDs{iName}, '_downstream_partners.csv'];
    
    % Load upstream partners
    us_partners = readtable(fullfile(parentDir, us_file_name), 'delimiter', ',');
    ID_table = repmat(table({curr_name}, neuron_IDs(iName), {'upstream'}, 'variableNames', ...
            {'neuron_name', 'neuron_ID', 'partnerDirection'}), size(us_partners, 1), 1);
    us_partners = [ID_table, us_partners];
        
    % Load downstream partners
    ds_partners = readtable(fullfile(parentDir, ds_file_name), 'delimiter', ',');
    ID_table = repmat(table({curr_name}, neuron_IDs(iName), {'downstream'}, 'variableNames', ...
            {'neuron_name', 'neuron_ID', 'partnerDirection'}), size(ds_partners, 1), 1);
    ds_partners = [ID_table, ds_partners];
    
    % Add both to main table
    tbl = [tbl; us_partners];
    tbl = [tbl; ds_partners];
end
    
% save_file = fullfile(parentDir, 'full_connectivity_table_backup.csv');
% writetable(tbl, save_file);

catch ME; rethrow(ME); end

%% LOAD SAVED CONNECTIVTY DATA
parentDir = 'C:\Users\Wilson Lab\Google Drive\Lab Work\LH-DAN_connectomics_analysis';
tbl = readtable(fullfile(parentDir, 'full_connectivity_table_all_cells.csv'), ...
        'delimiter', ',');

% DNs = readtable(fullfile(parentDir, 'putative_DN_IDs.csv'), 'delimiter', ',');    

%% Create tables of reciprocal connectivity

topPartnerCount = [];
saveTables = 0;


tracedTbl = tbl(strcmpi(tbl.partnerStatus, 'traced') & tbl.partnerSize > 1.2e8, :);
tracedTbl = tracedTbl(~strcmp(tracedTbl.partnerType, 'NA'), :);     

% Grouped by cell type
groupTbl = groupsummary(tracedTbl, {'neuronID', 'neuronName', 'neuronType', ...
        'partnerType', 'partnerDirection'}, 'sum', {'synapseCount'});
groupTbl.Properties.VariableNames{strcmp(groupTbl.Properties.VariableNames, ...
        'sum_synapseCount')} = 'synapseCount';
groupTbl.Properties.VariableNames{strcmp(groupTbl.Properties.VariableNames, ...
        'GroupCount')} = 'nCells';
tblDS = groupTbl(strcmp(groupTbl.partnerDirection, 'downstream'), :);
tblUS = groupTbl(strcmp(groupTbl.partnerDirection, 'upstream'), :);
mergeTbl = outerjoin(tblDS, tblUS, 'keys', {'neuronID', 'neuronName', 'neuronType', 'partnerType'}, ...
        'mergekeys', 1);
mergeTbl.synapseCount_tblUS(isnan(mergeTbl.synapseCount_tblUS)) = 0;
mergeTbl.nCells_tblUS(isnan(mergeTbl.nCells_tblUS)) = 0;
mergeTbl.synapseCount_tblDS(isnan(mergeTbl.synapseCount_tblDS)) = 0;
mergeTbl.nCells_tblDS(isnan(mergeTbl.nCells_tblDS)) = 0;
cellTypeRecipTable = mergeTbl(~isnan(mergeTbl.nCells_tblDS) & ~isnan(mergeTbl.nCells_tblUS), ...
        {'neuronID', 'neuronType', 'partnerType', 'nCells_tblDS', 'synapseCount_tblDS', ...
        'nCells_tblUS', 'synapseCount_tblUS'});
cellTypeRecipTable.totalSynapses = cellTypeRecipTable.synapseCount_tblDS + ...
        cellTypeRecipTable.synapseCount_tblUS;
if ~isempty(topPartnerCount)
    newTable = [];
    neuronList = unique(cellTypeRecipTable.neuronID);
    for iType = 1:numel(neuronList)
        currTypeTbl = cellTypeRecipTable(strcmp(cellTypeRecipTable.neuronID, neuronList{iType}), :);
        newTable = [newTable; tail(sortrows(currTypeTbl(:, {'neuronID', 'neuronType', 'partnerType', ...
            'nCells_tblDS', 'synapseCount_tblDS', 'nCells_tblUS', 'synapseCount_tblUS', ...
            'totalSynapses'}), 'totalSynapses'), topPartnerCount)];
    end
    cellTypeRecipTable = newTable;
end

% Single cell connections
singleCellTbl = tracedTbl(:, {'neuronID', 'neuronName', 'neuronType', 'partnerID', 'partnerType', ...
        'partnerDirection', 'synapseCount'});
tblDS = singleCellTbl(strcmp(singleCellTbl.partnerDirection, 'downstream'), :);
tblUS = singleCellTbl(strcmp(singleCellTbl.partnerDirection, 'upstream'), :);
mergeTbl = outerjoin(tblDS, tblUS, 'keys', {'neuronID', 'neuronName', 'neuronType', 'partnerID', ...
        'partnerType'}, 'mergekeys', 1);
singleCellRecipTable = mergeTbl(~isnan(mergeTbl.synapseCount_tblUS) & ...
        ~isnan(mergeTbl.synapseCount_tblDS), {'neuronID', 'neuronType', 'partnerID', 'partnerType', ...
        'synapseCount_tblDS', 'synapseCount_tblUS'});
singleCellRecipTable.totalSynapses = singleCellRecipTable.synapseCount_tblDS + ...
        singleCellRecipTable.synapseCount_tblUS;
if ~isempty(topPartnerCount)
    newTable = [];
    neuronList = unique(singleCellRecipTable.neuronID);
    for iType = 1:numel(neuronList)
        currTypeTbl = singleCellRecipTable(strcmp(singleCellRecipTable.neuronID, ...
                neuronList{iType}), :);
        newTable = [newTable; tail(sortrows(currTypeTbl(:, {'neuronID', 'neuronType', 'partnerID', ...
            'partnerType', 'synapseCount_tblDS', 'synapseCount_tblUS', 'totalSynapses'}), ...
            'totalSynapses'), topPartnerCount)];
    end
    singleCellRecipTable = newTable;
end

% Save tables
if saveTables
    writetable(cellTypeRecipTable, fullfile(parentDir, 'cellTypeReciprocityTable.csv'));
    writetable(singleCellRecipTable, fullfile(parentDir, 'singleCellReciprocityTable.csv'));
end


%
plotTable = cellTypeRecipTable(cellTypeRecipTable.totalSynapses > 1, :);
% plotTable = singleCellRecipTable;
plotTable.eqIndex = abs(0.5 - (plotTable.synapseCount_tblDS ./ ...
        plotTable.totalSynapses));
plotTable = sortrows(plotTable, 'totalSynapses', 'descend');
f = figure(2);clf;hold on;
% neuronList = sortrows(unique(cellTypeRecipTable(:, 1:2)), 2);
neuronList = sortrows(unique(plotTable(:, 2)), 1);
% cm = jet(size(neuronList, 1)) * 0.95;
cm = [rgb('blue'); rgb('blue'); rgb('orange'); rgb('red'); rgb('red'); rgb('red'); rgb('magenta'); ...
    rgb('magenta'); rgb('magenta');rgb('green'); rgb('violet');rgb('cyan'); rgb('cyan');rgb('cyan'); ...
    rgb('indigo');rgb('indigo'); rgb('indigo');rgb('indigo'); rgb('indigo');];
meanVals = [];
plotDims = numSubplots(size(neuronList,1));
for iCell = 1:size(neuronList, 1)
       ax = subaxis(plotDims(1), plotDims(2), iCell); cla();
%      yy = sort(plotTable.eqIndex(strcmp(cellTypeRecipTable.neuronID, ...
%             neuronList.neuronID{iCell})));
%      yy = (plotTable.eqIndex(strcmp(plotTable.neuronType, ...
%             neuronList.neuronType{iCell})));
%      xx = (1:numel(yy)) ./ numel(yy);
%      meanVals(iCell) = mean(yy');
     xx = plotTable.synapseCount_tblDS(strcmp(plotTable.neuronType, ...
            neuronList.neuronType{iCell}));
     yy = plotTable.synapseCount_tblUS(strcmp(plotTable.neuronType, ...
            neuronList.neuronType{iCell}));
%      disp([neuronList.neuronType{iCell}, ': ', num2str(mean(yy'), 3)])
     plot(ax, xx, yy, 'x', 'color', cm(iCell, :));   
     hold on;
     xL = xlim();
     yL = ylim();
     plot(ax, [1, max([xL, yL])], [1, max([xL, yL])], '-', 'color', 'k')
     title(regexprep(neuronList.neuronType{iCell}, '_', '\\_'));
     xlabel('Output synapses');
     ylabel('Input synapses');
     ax.XScale = 'log'; 
     ax.YScale = 'log';
end
% neuronList.meanVals = meanVals';
% disp(sortrows(neuronList, 2, 'descend'))
% legend(neuronList.neuronType, 'location', 'nw');

%% 
dt = DataTable(tbl);
filterDefs = struct();
typeList = unique(dt.sourceData.neuronType);

currType = 8;

% filterDefs.DAN_name = 'TypeB';
filterDefs.neuronType = typeList{currType};
filterDefs.partnerDirection = 'upstream';

thresh = 5;

% Get data subset
dt = dt.initialize_filters(filterDefs);
sub = dt.subset;

% Drop rows for partners that are incompletely traced or that do not have an assigned cell type
sub = sub(strcmpi(sub.status, 'traced') & ~cellfun(@isempty, sub.partnerType) & ...
        ~strcmp(sub.partnerType, 'NA'), :);

% Calculate connection weight of with each partner neuron as a percentage of total DAN synapses
totalSynapses = sum(sub.synapseCount);
connectionWeightPct = 100 * (sub.synapseCount ./ totalSynapses);
sub.connectionWeightPct = connectionWeightPct;

% Create another table that groups partners by cell type
groupedSub = groupsummary(sub, {'neuronID', 'neuronName', 'partnerDirection', 'partnerType'}, ...
        'sum', {'synapseCount'});
groupedSub.connectionWeightPct = 100 * (groupedSub.sum_synapseCount ./ totalSynapses);
    
% Plot figures
f = figure(1);clf;
f.Color = [1 1 1];
colororder(f, [0, 0, 0; 0, 0, 0]);
if strcmp(filterDefs.partnerDirection, 'downstream')
    yLabelStr = 'output';
else
    yLabelStr = 'input';
end

% First plot percentages grouped by cell types
ax = subaxis(2, 1, 1, 'ml', 0.14, 'mr', 0.14, 'sv', 0.07, 'mb', 0.07, 'mt', 0.04);
xx = (1:numel(sort(groupedSub.connectionWeightPct))) ./ numel(sort(groupedSub.connectionWeightPct)); 
plot(cumsum(sort(groupedSub.connectionWeightPct)), '-', 'linewidth', 1.5);
ax.Children.Marker = '.';
ax.Children.MarkerSize = 12;
ax.FontSize = 13;

% Set axes limits on a log scale
maxVal = max(groupedSub.connectionWeightPct);
minVal = min(groupedSub.connectionWeightPct);
ax = logify(ax, [log10(exp_floor(minVal)), log10(exp_ceil(maxVal))]);
ax.YLim = [minVal, maxVal] .* [0.9, 1.5];
newYL = ax.YLim;
ylabel(['% of total ', yLabelStr, ' synapses']);

% usMaxVals(currType) = maxVal;


% Add a matching raw synapse count scale to other side of the plot 
yyaxis(ax, 'right')
ax.FontSize = 13;
ylabel('Synapse count');
title([regexprep(filterDefs.neuronType, '_', '\\_'), ' ', filterDefs.partnerDirection, ...
        ' partner cell types'], 'fontsize', 14)
ax.YLim = (newYL / 100) .* totalSynapses; % Match left and right axes limits
logify(ax, [log10(exp_floor(ax.YLim(1))), log10(exp_ceil(ax.YLim(2)))]);

% Now add another version of the same plot using single neurons instead of cell types
ax = subaxis(2, 1, 2);
plot(sort(sub.connectionWeightPct), '-', 'linewidth', 1.5);
ax.Children.Marker = '.';
ax.Children.MarkerSize = 12;
ax.FontSize = 13;

% Set axes limits on a log scale
maxVal = max(sub.connectionWeightPct);
minVal = min(sub.connectionWeightPct);
ax = logify(ax, [log10(exp_floor(minVal)), log10(exp_ceil(maxVal))]);
ax.YLim = [minVal, maxVal] .* [0.9, 1.5];
newYL = ax.YLim;
ylabel(['% of total ', yLabelStr, ' synapses']);
xlabel('Cell or cell type number')
title(['Individual ', filterDefs.partnerDirection, ' partner neurons'], 'fontsize', 14)

% Add a matching raw synapse count scale to other side of the plot 
yyaxis(ax, 'right')
ax.FontSize = 13;
ylabel('Synapse count');
ax.YLim = (newYL / 100) .* totalSynapses; % Match left and right axes limits
ax = logify(ax, [log10(exp_floor(ax.YLim(1))), log10(exp_ceil(ax.YLim(2)))]);

% Create summary tables of the top partner cells and cell types
disp(' '); disp(' ')

nRows = 20;

cellTypeTbl = tail(sortrows(groupedSub, 'connectionWeightPct'), nRows);
cellTypeTbl.Properties.VariableNames{strcmp(cellTypeTbl.Properties.VariableNames, ...
        'sum_synapseCount')} = 'synapseCount';
cellTypeTbl.Properties.VariableNames{strcmp(cellTypeTbl.Properties.VariableNames, ...
        'GroupCount')} = 'nCells';
    
singleCellTbl = tail(sortrows(sub(:, {'neuronID', 'neuronName', 'partnerName', 'partnerType', ...
        'partnerDirection', 'synapseCount', 'connectionWeightPct'}), 'connectionWeightPct'), nRows);
    
disp(cellTypeTbl)
disp(singleCellTbl);

%% Local functions 

% Change Y axis to a log scale and set tick labels
function ax = logify(ax, yExpLims)

YTickVals = [];
baseYTickVals = [10^yExpLims(1):10^yExpLims(1):(10^yExpLims(1) * 5)];
for i = 1:diff(yExpLims)
    YTickVals = [YTickVals, baseYTickVals * 10^(i-1)];
end
YTickVals(end + 1) = 10^yExpLims(2);
ax.YScale = 'log';
ax.YTick = YTickVals;

end

% Round up to the next power of 10
function x = exp_ceil(n, startingExp)
    if nargin < 2
        startingExp = -5;
    end
    if n > 10^startingExp
        x = exp_ceil(n, startingExp + 1);        
    else
        x = 10^startingExp;
    end
end

% Round down to the next power of 10
function x = exp_floor(n, startingExp)
    if nargin < 2
        startingExp = 10;
    end
    if n < 10^startingExp
        x = exp_floor(n, startingExp - 1);        
    else
        x = 10^startingExp;
    end
end

% function

%%
% f = figure(2);clf;
% f.Color = [1 1 1];
% ax = subaxis(1, 2, 1);
% plot(sort(groupedSub.sum_synapseCount), '-');
% title('Cell types', 'fontsize', 14)
% ax.YScale = 'log';
% ax.YTick = [1:5, 20:10:50, 100, 200:100:500, 1000];
% ax = subaxis(1, 2, 2);
% plot(sort(sub.synapseCount), '-');
% title('Single neurons', 'fontsize', 14)
% ax.YScale = 'log';
% ax.YTick = [1:5, 20:10:50, 100, 200:100:500, 1000];
% suptitle('Connection synapse counts')

% figure(3);clf; hold on;
% plot(sort(groupedSub.partnerConnectionWeightPct), '-o');
% plot(sort(sub.partnerConnectionWeightPct), '-o')
% legend({'Cell types', 'Single neurons'}, 'fontsize', 14, 'location', 'nw')
% title('Connection weights (% of total partner neuron input/outputs)')
% 
% Show tables with top partners

% 
% % Display percentage of partners with connection weight above a certain threshold
% nStrongPartners = sum(groupedSub.sum_synapseCount > thresh);
% pctStrongPartners = 100 * (nStrongPartners / numel(groupedSub.sum_synapseCount));
% disp([num2str(round(pctStrongPartners, 1)), '% of ', filterDefs.partnerDirection, ...
%         ' cell types connected by >', num2str(thresh), ' synapses']);
% nStrongPartners = sum(sub.synapseCount > thresh);
% pctStrongPartners = 100 * (nStrongPartners / numel(sub.synapseCount));
% disp([num2str(round(pctStrongPartners, 1)), '% of ', filterDefs.partnerDirection, ...
%         ' neurons connected by >', num2str(thresh), ' synapses']);


% figure(4);clf; plot(groupedSub.sum_partnerOutputSynapseCount, groupedSub.sum_synapseCount, '.')
% xlabel('partner synapses')
% ylabel('connection strength (synapses)')



























    

