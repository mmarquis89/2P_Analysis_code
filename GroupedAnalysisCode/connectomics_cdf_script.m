dt = DataTable(tbl);
filterDefs = struct();

% Combining data for the two halves of the PPL1 DAN with the cropped contralateral projection
dt.sourceData(strcmp(dt.sourceData.neuronID, '642516892'), 'neuronID') = {'5813019513'};
dt.sourceData(strcmp(dt.sourceData.neuronName, 'PPL105(a''2a2)(PDL05)_L'), 'neuronName') = {'PPL105(a''2a2)_R'};

typeList = sortrows(unique(dt.sourceData(:, {'neuronID', 'neuronName'})), 2);
% typeList = unique(dt.sourceData.neuronType);
neuronList = unique(dt.sourceData.neuronID);
neuronNameList = unique(dt.sourceData.neuronName);

currTypeList = typeList([1:13, 19:21, 15:18],:);
% currTypeList = typeList([1:15, 20:22, 16:19], :);
partnerDirection = 'upstream';
yVarType =  'percentage'; %'count'; %
yScaleType =  'linear'; %'linear'; % 
xNorm = 1;
xMinPct = 80;
groupPartnersByType = 1;
tracedPartnersOnly = 1;
legendPos = 'nw';
saveFig = 0;

% cm = jet(size((currTypeList), 1)) * 0.95;

cm = [rgb('blue'); rgb('blue'); ...
        rgb('cyan'); rgb('cyan'); rgb('cyan'); ...
        rgb('red'); rgb('red'); rgb('red'); ...
        rgb('green'); rgb('green'); rgb('green'); ...
        rgb('magenta'); rgb('magenta'); ...
%         rgb('black'); ...
        rgb('orange'); rgb('orange'); rgb('orange'); ...
        rgb('brown'); rgb('brown'); rgb('brown'); rgb('brown')];


% Create figure
f = figure(1);clf;
f.Color = [1 1 1];
ax = axes(); hold on;

filterDefs.partnerDirection = partnerDirection;
if strcmp(filterDefs.partnerDirection, 'downstream')
    directionStr = 'output';
else
    directionStr = 'input';
end

partnerTypeCounts = [];
maxSynapseCounts = [];
totalSynapseCounts = [];

overallMinVal = inf;
overallMaxVal = 0;
cellTypeTbl_allTypes = [];
singleCellTbl_allTypes = [];
count = 1;
for iType = 1:size((currTypeList), 1)
    
%     filterDefs.neuronType = currTypeList{iType};
    filterDefs.neuronID = currTypeList.neuronID{iType};

    % Get data subset
    dt = dt.initialize_filters(filterDefs);
    sub = dt.subset;
    
    % Drop rows for partners that are incompletely traced or too small
    if tracedPartnersOnly
        sub = sub(strcmpi(sub.partnerStatus, 'traced') & sub.partnerSize > 1.2e8, :);
    end
        
    % Calculate connection weight with each partner neuron as a percentage of total synapses
    totalSynapses = sum(sub.synapseCount);
    maxTotalSynapses(iType) = max(totalSynapses);
    connectionWeightPct = 100 * (sub.synapseCount ./ totalSynapses);
    sub.connectionWeightPct = connectionWeightPct;
    
    % Assign a placeholder cell type to each cell that doesn't have one
    for iPartner = 1:size(sub, 1)
        if strcmp(sub.partnerType{iPartner}, 'NA')
            sub.partnerType{iPartner} = ['unknown_', pad(num2str(count), 4, 'left', '0')];
            count = count + 1;
        end
    end
    
    % Create another table that groups partners by cell type
    groupedSub = groupsummary(sub, {'neuronID', 'neuronName', 'partnerDirection', 'partnerID', ...
            'partnerType'}, ...
            'sum', {'synapseCount'});
%     groupedSub = groupedSub(~strcmp(groupedSub.partnerType, 'NA'), :);
    groupedSub.connectionWeightPct = 100 * (groupedSub.sum_synapseCount ./ totalSynapses);
    
    totalSynapseCounts(iType) = sum(groupedSub.sum_synapseCount);
    partnerTypeCounts(iType) = numel(groupedSub.sum_synapseCount);
    maxSynapseCounts(iType) = max(groupedSub.sum_synapseCount);
    
    % Get correct plot data for selected options
    if groupPartnersByType
        pctData = groupedSub.connectionWeightPct;
        countData = groupedSub.sum_synapseCount;
    else
        pctData = sub.connectionWeightPct;
        countData = sub.synapseCount;
    end
    if strcmpi(yVarType, 'count')
        yy = sort(countData);
    else
        yy = sort(pctData);
    end
    if xNorm
        xx = (1:numel(yy)) ./ numel(yy);
    else
        xx = 1:numel(yy);
    end
    
%     xx = xx - (numel(xx) - 100);
    
    % Keep track of overall min and max vals for setting final axes limits
    overallMinVal = min(overallMinVal, min(yy));
    overallMaxVal = max(overallMaxVal, max(yy));
    
%     yy = smoothdata(yy, 'gaussian', 3);
%     yy = cumsum(yy);
    
    % Plot data
    plot(ax, xx, yy, '-', 'linewidth', 1.5, 'color', cm(iType, :));
    ax.Children(1).Marker = '.';
    ax.Children(1).MarkerSize = 12;
    ax.FontSize = 14; 
    if iType == 1
        lineObjects = ax.Children;
    else
        lineObjects(end + 1) = ax.Children(1);
    end
    
    % Create summary tables of the top partner cells and cell types    
    nRows = 15;
    
    cellTypeTbl = tail(sortrows(groupedSub, 'connectionWeightPct'), nRows);
    cellTypeTbl.Properties.VariableNames{strcmp(cellTypeTbl.Properties.VariableNames, ...
            'sum_synapseCount')} = 'synapseCount';
    cellTypeTbl.Properties.VariableNames{strcmp(cellTypeTbl.Properties.VariableNames, ...
            'GroupCount')} = 'nCells';
    
    singleCellTbl = tail(sortrows(sub(:, {'neuronID', 'neuronName', 'partnerName', 'partnerType', ...
            'partnerDirection', 'synapseCount', 'connectionWeightPct'}), 'connectionWeightPct'), ...
            nRows);
    
    cellTypeTbl_allTypes = [cellTypeTbl_allTypes; cellTypeTbl];
    singleCellTbl_allTypes = [singleCellTbl_allTypes; singleCellTbl];
    
end

% Change Y-axis to a log scale if necessary
if strcmp(yScaleType, 'log')
    ax = logify(ax, [log10(exp_floor(overallMinVal)), log10(exp_ceil(overallMaxVal))]);
end

% Set Y limits
if strcmp(yScaleType, 'log')
    ax.YLim = [overallMinVal, overallMaxVal] .* [0.9, 1.5];
else
    ax.YLim = [overallMinVal, overallMaxVal] .* [0.9, 1.1];
end

% Set X limits if necessary
if xNorm
    xlim(ax, [xMinPct/100, 1.02]);
%     ax.XTickLabels = '';
elseif xMinPct > 0
    xL = xlim();
    pctThresh = round(xL(2) * (xMinPct / 100));
    xlim([pctThresh, xL(2)])
end

% Add legend and labels
if strcmp(yVarType, 'count')
    ylabel('Synapse count');
else
    ylabel(['% of total ', directionStr, ' synapses']);
end
if groupPartnersByType
    xlabel([capitalize(partnerDirection), ' partner cell types'])
else
    xlabel([capitalize(partnerDirection), ' partner cells'])
end
if xNorm
    ax.XLabel.String = [ax.XLabel.String, ' (percentile)'];
end
legend(lineObjects, regexprep(currTypeList.neuronName, '_', '\\_'), 'location', legendPos)
title([capitalize(partnerDirection), ' synaptic partner connection strengths'])

if saveFig
    saveDir = fullfile(parentDir, 'figs');
    
    if groupPartnersByType
        partnerTypeStr = 'cellTypes';
    else
        partnerTypeStr = 'cells';
    end
    if xNorm
        xNormStr = 'normalized';
    else
        xNormStr = 'raw';
    end
    
    fileName = [partnerDirection, '_partner_', partnerTypeStr, '_synapse', capitalize(yVarType), ...
        's_', yScaleType, 'Y_', xNormStr, 'X_xMin_', num2str(xMinPct)];
    
    save_figure(f, saveDir, fileName);
end

% % Save tables with top partner cells and cell types
disp(cellTypeTbl_allTypes);
% disp(singleCellTbl_allTypes);
% 
% fileName = ['top_', partnerDirection, '_partner_cellTypes.csv'];
% writetable(cellTypeTbl_allTypes, fullfile(parentDir, fileName));
% 
% fileName = ['top_', partnerDirection, '_partner_cells.csv'];
% writetable(singleCellTbl_allTypes, fullfile(parentDir, fileName));
% 
%%
% % tblDS = cellTypeTbl_allTypes;
% 
% test = outerjoin(tblUS, tblDS, 'Keys', {'neuronID', 'partnerID'}, 'MergeKeys', 1, 'LeftVariables', ...
%         {'neuronID', 'neuronName', 'partnerID', 'partnerType', 'synapseCount'}, 'RightVariables', ...
%         {'neuronID', 'neuronName', 'partnerID', 'partnerType', 'synapseCount'});
% 
% test2 = test(~isnan(test.synapseCount_tblUS) & ~isnan(test.synapseCount_tblDS), :);
% 
% writetable(test2, fullfile('C:\Users\Wilson Lab\Google Drive\Lab Work\LH-DAN_connectomics_analysis', ...
%         'testRecipTable.csv'));
%% Create scatter plots of summary variables
tb = table(currTypeList.neuronName, totalSynapseCounts', partnerTypeCounts', maxSynapseCounts', ...
        'variablenames', {'cellType', 'totalSynapses', 'totalPartnerTypes', 'maxSynapseCount'});

xx = tb.totalSynapses;
xLab = ['Total ', directionStr, ' synapses'];

yy = tb.maxSynapseCount;
yLab = ['Max synapses with a single partner cell type'];
fileName = ['total_', directionStr, '_synapses_vs_max_synapses_with_single_cell_type'];
% 
yy = tb.totalPartnerTypes;
yLab = ['Number of ', partnerDirection, ' cell types'];
fileName = ['total_', directionStr, '_synapses_vs_total_', partnerDirection, '_cell_types'];

saveFig = 0;

dx = ones(size(xx)) * 100;
% dx(7) = dx(7) - 850;
dy = ones(size(xx)) * 25;
% dy(9) = dy(9) - 50;

f = figure(4);clf;
f.Color = [1 1 1];
hold on; 
plot(xx, yy, '.', 'color', 'b', 'markersize', 25);  
ax = gca();
ax.FontSize = 14;
xlabel(xLab)
ylabel(yLab)
c = regexprep(tb.cellType, '_', '\\_');
text(xx + dx, yy + dy, c, 'FontSize', 11);

% labelpoints(xx([1:2, 4:13, 16]), yy([1:2, 4:13, 16]), c([1:2, 4:13, 16]), 'fontsize', 11, 'adjust_axes', 1, 'buffer', 0.1, 'position', 'nw')
% labelpoints(xx([3, 14]), yy([3, 14]), c([3, 14]), 'fontsize', 11, 'adjust_axes', 1, 'buffer', 0.1, 'position', 'ne')
% labelpoints(xx([15]), yy([15]), c([15]), 'fontsize', 11, 'adjust_axes', 1, 'buffer', 0.1, 'position', 'sw')

if saveFig
    save_figure(f, saveDir, fileName);
end

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
